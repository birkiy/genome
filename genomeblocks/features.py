"""Feature container types: Features, Tags."""
from typing import Any

class Features(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __str__(s):
        lines = [f"{s.__class__.__name__} (n={len(s)})"]
        lines.append(s.table())
        return "\n".join(lines)
    __repr__ = __str__

    def table(self):
        return "(No table representation for base Features)"


class Tags(Features):
    """a dictionary of tags"""
    def table(self):
        lines = []
        if len(self) > 0:
            lines.append(f"{'Name':<50} {'Count':<10}")
            lines.append("-" * 60)
            for k in self:
                lines.append(f"{k:<50} {len(self[k]):<10}")
        return "\n".join(lines)
