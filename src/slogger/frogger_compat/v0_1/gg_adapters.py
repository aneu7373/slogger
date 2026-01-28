from typing import List, Dict


def add_adapters(fragments: List[Dict], gg_cfg: List[Dict]) -> List[Dict]:
    """
    gg_cfg: list of adapter dicts indexed by frag_index-1
      { left_adapter: "...", right_adapter: "..." }
    """
    out = []
    for f in fragments:
        idx = int(f["frag_index"]) - 1
        if idx < 0 or idx >= len(gg_cfg):
            raise ValueError(f"No Golden Gate adapter specified for frag_index={f['frag_index']}")
        left = gg_cfg[idx]["left_adapter"]
        right = gg_cfg[idx]["right_adapter"]
        out.append({**f, "cloning_seq": f"{left}{f['core_seq']}{right}"})
    return out
