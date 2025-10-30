"""Tag management utilities for annotating `Loci` collections."""
from __future__ import annotations

import operator
import re
from dataclasses import dataclass
from typing import Any, Callable, Dict, Iterable, Mapping, Optional, Sequence, Set, Union, TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from .loci import Loci

_NON_ALPHA = re.compile(r"[^0-9a-zA-Z]+")


@dataclass
class TagEntry:
    """Internal representation of a stored tag."""

    key: str
    kind: Literal["set", "mapping"]
    data: Union[frozenset[str], Dict[str, Any]]


class Tags:
    """In-memory tag store with expression-based querying."""

    def __init__(self, loci: "Loci", *, verbose: bool = True):
        from .loci import Loci  # Delayed import to avoid eager circular dependency

        if not isinstance(loci, Loci):
            loci = Loci(loci)
        self._loci = loci
        self._uid_index = {l.uid: idx for idx, l in enumerate(self._loci)}
        self._tags: Dict[str, TagEntry] = {}
        self._alias_index: Dict[str, Set[str]] = {}
        self._key_aliases: Dict[str, Set[str]] = {}
        self._verbose = verbose

    @classmethod
    def make(cls, loci: "Loci", *, verbose: bool = True) -> "Tags":
        """Construct a tag container anchored to the provided loci."""
        tags = cls(loci, verbose=verbose)
        tags._log(f"[INFO] Seeded Tags with {len(tags._loci)} loci.")
        return tags

    def add(self, mapping: Mapping[str, Any]) -> "Tags":
        """Register tag definitions from a mapping of name -> Loci/UIDs/dict."""
        if not isinstance(mapping, Mapping):
            raise TypeError("Tags.add expects a mapping of tag names to payloads.")

        for key, payload in mapping.items():
            entry = self._build_entry(key, payload)
            self._register_aliases(entry)
            self._tags[key] = entry
            self._log(f"[INFO] Stored tag '{key}' ({entry.kind}).")
        return self

    def query(self, predicate: Callable[["TagAccessor"], "TagSelection"]) -> "Loci":
        """Evaluate an expression over tags and return the matching loci."""
        accessor = TagAccessor(self)
        result = predicate(accessor)
        selection = self._coerce_query_result(result)
        return self._subset_from_uids(selection.uids)

    def __contains__(self, key: str) -> bool:  # pragma: no cover - trivial
        try:
            self._locate_key(key)
        except AttributeError:
            return False
        return True

    def __getitem__(self, key: str) -> Union["TagSelection", "TagNumericView"]:
        resolved = self._locate_key(key)
        return self._view_for_key(resolved)

    def keys(self) -> Iterable[str]:  # pragma: no cover - trivial
        return self._tags.keys()

    def items(self) -> Iterable[tuple[str, Union[frozenset[str], Dict[str, Any]]]]:
        for key, entry in self._tags.items():
            yield key, entry.data

    def uids(self) -> Sequence[str]:
        """Return the ordered list of reference locus UIDs."""
        return list(self._uid_index.keys())

    def table(self) -> str:
        """Simple textual summary of stored tags."""
        if not self._tags:
            return ""
        lines = [f"{'Name':<30} {'Type':<10} {'Count':<10}"]
        lines.append("-" * 52)
        for key in sorted(self._tags):
            entry = self._tags[key]
            count = len(entry.data)
            lines.append(f"{key:<30} {entry.kind:<10} {count:<10}")
        return "\n".join(lines)
    
    def __str__(self):
        return self.table()
    __repr__ = __str__

    # Internal helpers -------------------------------------------------

    def _log(self, message: str) -> None:
        if self._verbose:
            print(message)

    def _build_entry(self, key: str, payload: Any) -> TagEntry:
        if isinstance(payload, Mapping):
            data = self._prepare_mapping(payload, key)
            return TagEntry(key=key, kind="mapping", data=data)
        if self._is_loci_instance(payload):
            uids = self._uids_from_loci(payload)
            return TagEntry(key=key, kind="set", data=frozenset(uids))
        if self._is_uid_iterable(payload):
            uids = self._validate_uids((str(uid) for uid in payload), key)
            return TagEntry(key=key, kind="set", data=frozenset(uids))
        raise TypeError(
            f"Unsupported payload for tag '{key}'. Expected Loci, iterable of UIDs, or Mapping.",
        )

    def _prepare_mapping(self, mapping: Mapping[Any, Any], key: str) -> Dict[str, Any]:
        cleaned: Dict[str, Any] = {}
        for raw_uid, raw_value in mapping.items():
            uid = str(raw_uid)
            if uid not in self._uid_index:
                raise KeyError(f"Tag '{key}' references unknown uid '{uid}'.")
            cleaned[uid] = self._coerce_mapping_value(raw_value)
        return cleaned

    @staticmethod
    def _coerce_mapping_value(value: Any) -> Any:
        if isinstance(value, (int, float)):
            return value
        if isinstance(value, str):
            try:
                return float(value)
            except ValueError:
                return value
        return value

    def _register_aliases(self, entry: TagEntry) -> None:
        if entry.key in self._key_aliases:
            for alias in self._key_aliases[entry.key]:
                holders = self._alias_index.get(alias)
                if holders:
                    holders.discard(entry.key)
                    if not holders:
                        self._alias_index.pop(alias)
        aliases = self._make_aliases(entry.key)
        self._key_aliases[entry.key] = aliases
        for alias in aliases:
            self._alias_index.setdefault(alias, set()).add(entry.key)

    @staticmethod
    def _make_aliases(name: str) -> Set[str]:
        lowered = name.lower()
        aliases = {lowered}
        tokens = [token for token in _NON_ALPHA.split(lowered) if token]
        aliases.update(tokens)
        compact = _NON_ALPHA.sub("_", lowered).strip("_")
        if compact:
            aliases.add(compact)
        joined = "".join(tokens)
        if joined:
            aliases.add(joined)
        return aliases

    def _locate_key(self, alias: str) -> str:
        if alias in self._tags:
            return alias
        lowered = alias.lower()
        for key in self._tags:
            if key.lower() == lowered:
                return key
        matches = self._alias_index.get(lowered)
        if matches:
            if len(matches) == 1:
                return next(iter(matches))
            raise AttributeError(
                f"Ambiguous tag alias '{alias}' matches {sorted(matches)}. Use the full tag name.",
            )
        raise AttributeError(f"Unknown tag '{alias}'. Available tags: {sorted(self._tags)}")

    def _view_for_key(self, key: str) -> Union["TagSelection", "TagNumericView"]:
        entry = self._tags[key]
        if entry.kind == "set":
            return TagSelection(self, entry.data, name=key)
        return TagNumericView(self, key, entry.data)

    def _coerce_query_result(
        self,
        value: Any,
    ) -> "TagSelection":
        if isinstance(value, PendingNumericCombination):
            raise ValueError("Numeric tag comparison is missing a comparator.")
        if isinstance(value, TagNumericView):
            raise ValueError("Numeric tags must be compared to a value (e.g. l.atac > 1.5).")
        selection = self._maybe_selection(value)
        if selection is None:
            raise TypeError("Query must return a tag selection or iterable of known UIDs.")
        return selection

    def _maybe_selection(self, value: Any) -> Optional["TagSelection"]:
        if isinstance(value, TagSelection):
            return value
        if isinstance(value, PendingNumericCombination):
            return None
        if isinstance(value, TagNumericView):
            return None
        if self._is_loci_instance(value):
            uids = self._uids_from_loci(value)
            return TagSelection(self, uids)
        if self._is_uid_iterable(value):
            uids = self._validate_uids((str(uid) for uid in value))
            return TagSelection(self, uids)
        return None

    def _subset_from_uids(self, uids: Iterable[str]) -> "Loci":
        keep = set(uids)
        from .loci import Loci  # Delayed to avoid top-level cycle

        return Loci([locus for locus in self._loci if locus.uid in keep])

    def _uids_from_loci(self, loci: Any) -> Set[str]:
        from .loci import Loci  # Delayed import

        if not isinstance(loci, Loci):
            loci = Loci(loci)
        return self._validate_uids((l.uid for l in loci))

    def _validate_uids(self, uids: Iterable[str], key: Optional[str] = None) -> Set[str]:
        resolved: Set[str] = set()
        for uid in uids:
            if uid not in self._uid_index:
                if key:
                    raise KeyError(f"Tag '{key}' references unknown uid '{uid}'.")
                raise KeyError(f"Unknown uid '{uid}'.")
            resolved.add(uid)
        return resolved

    @staticmethod
    def _is_uid_iterable(value: Any) -> bool:
        if isinstance(value, (str, bytes)):
            return False
        return isinstance(value, Iterable)

    @staticmethod
    def _is_loci_instance(value: Any) -> bool:
        try:
            from .loci import Loci  # Local import to avoid circular dependency
        except Exception:  # pragma: no cover - defensive
            return False
        return isinstance(value, Loci)


class TagAccessor:
    """Provides attribute-based access to tags inside query lambdas."""

    def __init__(self, tags: Tags):
        self._tags = tags

    def __getattr__(self, name: str) -> Union["TagSelection", "TagNumericView"]:
        key = self._tags._locate_key(name)
        return self._tags._view_for_key(key)

    def __getitem__(self, key: str) -> Union["TagSelection", "TagNumericView"]:
        return self.__getattr__(key)


class TagSelection:
    """Represents a set of locus UIDs backed by the parent Tags instance."""

    def __init__(self, tags: Tags, uids: Iterable[str], *, name: Optional[str] = None):
        self._tags = tags
        self._uids = frozenset(uids)
        self._name = name

    def __and__(self, other: Any) -> Union["TagSelection", "PendingNumericCombination"]:
        if isinstance(other, PendingNumericCombination):
            return other.__rand__(self)
        if isinstance(other, TagNumericView):
            return PendingNumericCombination(self, other)
        selection = self._tags._maybe_selection(other)
        if selection is None:
            raise TypeError("Cannot intersect TagSelection with unsupported type.")
        return TagSelection(self._tags, self._uids & selection._uids)

    def __rand__(self, other: Any) -> Union["TagSelection", "PendingNumericCombination"]:
        return self.__and__(other)

    def __or__(self, other: Any) -> "TagSelection":
        selection = self._tags._maybe_selection(other)
        if selection is None:
            raise TypeError("Cannot union TagSelection with unsupported type.")
        return TagSelection(self._tags, self._uids | selection._uids)

    def __ror__(self, other: Any) -> "TagSelection":
        return self.__or__(other)

    def __sub__(self, other: Any) -> "TagSelection":
        selection = self._tags._maybe_selection(other)
        if selection is None:
            raise TypeError("Cannot subtract unsupported type from TagSelection.")
        return TagSelection(self._tags, self._uids - selection._uids)

    def __rsub__(self, other: Any) -> "TagSelection":
        selection = self._tags._maybe_selection(other)
        if selection is None:
            raise TypeError("Cannot subtract TagSelection from unsupported type.")
        return TagSelection(self._tags, selection._uids - self._uids)

    def __xor__(self, other: Any) -> "TagSelection":
        selection = self._tags._maybe_selection(other)
        if selection is None:
            raise TypeError("Cannot xor TagSelection with unsupported type.")
        return TagSelection(self._tags, self._uids ^ selection._uids)

    def __rxor__(self, other: Any) -> "TagSelection":
        return self.__xor__(other)

    def __len__(self) -> int:  # pragma: no cover - trivial
        return len(self._uids)

    def __iter__(self):
        return iter(self._uids)

    def __contains__(self, uid: str) -> bool:  # pragma: no cover - trivial
        return uid in self._uids

    def __repr__(self) -> str:  # pragma: no cover - debug helper
        label = f" '{self._name}'" if self._name else ""
        return f"<TagSelection{label} n={len(self._uids)}>"

    @property
    def uids(self) -> Set[str]:
        return set(self._uids)

    def to_loci(self) -> "Loci":
        return self._tags._subset_from_uids(self._uids)


class TagNumericView:
    """Read-only view over a numeric/categorical tag mapping."""

    def __init__(self, tags: Tags, key: str, values: Mapping[str, Any]):
        self._tags = tags
        self._key = key
        self._values = dict(values)

    def __and__(self, other: Any) -> "PendingNumericCombination":
        selection = self._tags._maybe_selection(other)
        if selection is None:
            raise TypeError("Numeric tags can only be combined with selections using '&'.")
        return PendingNumericCombination(selection, self)

    def __rand__(self, other: Any) -> "PendingNumericCombination":
        return self.__and__(other)

    def __gt__(self, other: Any) -> TagSelection:
        return self._compare(operator.gt, other, ">")

    def __ge__(self, other: Any) -> TagSelection:
        return self._compare(operator.ge, other, ">=")

    def __lt__(self, other: Any) -> TagSelection:
        return self._compare(operator.lt, other, "<")

    def __le__(self, other: Any) -> TagSelection:
        return self._compare(operator.le, other, "<=")

    def __eq__(self, other: Any) -> TagSelection:  # type: ignore[override]
        return self._compare(operator.eq, other, "==")

    def __ne__(self, other: Any) -> TagSelection:  # type: ignore[override]
        return self._compare(operator.ne, other, "!=")

    def __getitem__(self, uid: str) -> Any:
        return self._values[uid]

    def keys(self) -> Iterable[str]:  # pragma: no cover - trivial
        return self._values.keys()

    def items(self) -> Iterable[tuple[str, Any]]:  # pragma: no cover - trivial
        return self._values.items()

    def _compare(self, op: Callable[[Any, Any], bool], other: Any, symbol: str) -> TagSelection:
        matches = {
            uid
            for uid, value in self._values.items()
            if self._safe_compare(op, value, other)
        }
        return TagSelection(self._tags, matches, name=f"{self._key} {symbol} {other}")

    @staticmethod
    def _safe_compare(op: Callable[[Any, Any], bool], left: Any, right: Any) -> bool:
        try:
            return bool(op(left, right))
        except TypeError:
            return False


class PendingNumericCombination:
    """Deferred intersection of a selection with a numeric tag comparison."""

    def __init__(self, base: TagSelection, numeric: TagNumericView):
        self._base = base
        self._numeric = numeric

    def __and__(self, other: Any) -> "PendingNumericCombination":
        selection = self._numeric._tags._maybe_selection(other)
        if selection is None:
            raise TypeError("Pending comparisons can only combine with TagSelections.")
        combined = self._base & selection
        if isinstance(combined, PendingNumericCombination):
            # Should not happen because base & selection returns TagSelection
            raise TypeError("Chained numeric comparisons without resolution are unsupported.")
        return PendingNumericCombination(combined, self._numeric)

    def __rand__(self, other: Any) -> "PendingNumericCombination":
        return self.__and__(other)

    def __gt__(self, other: Any) -> TagSelection:
        return self._resolve(operator.gt, other, ">")

    def __ge__(self, other: Any) -> TagSelection:
        return self._resolve(operator.ge, other, ">=")

    def __lt__(self, other: Any) -> TagSelection:
        return self._resolve(operator.lt, other, "<")

    def __le__(self, other: Any) -> TagSelection:
        return self._resolve(operator.le, other, "<=")

    def __eq__(self, other: Any) -> TagSelection:  # type: ignore[override]
        return self._resolve(operator.eq, other, "==")

    def __ne__(self, other: Any) -> TagSelection:  # type: ignore[override]
        return self._resolve(operator.ne, other, "!=")

    def _resolve(self, op: Callable[[Any, Any], bool], other: Any, symbol: str) -> TagSelection:
        selection = self._numeric._compare(op, other, symbol)
        return self._base & selection
