#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Aplica patches de compatibilidade pandas>=2 no PhaBox2 instalado (IDEMPOTENTE).

Por quê: PhaBox2 (>=2.0) foi escrito para idiomas do pandas<2. Com pandas>=2,
atribuicoes do tipo `df.loc[...] = "texto"` em colunas int/float levantam
`TypeError: Invalid value ... for dtype` / `LossySetitemError`, abortando
phagcn/cherry/votu. O env usa Python 3.12, onde pandas<2 nao tem wheel — por
isso corrigimos o codigo-fonte do PhaBox2 em vez de fazer downgrade do pandas.

Roda DENTRO do env conda do phabox (que tem phabox2 importavel). Seguro rodar
toda vez: se o patch ja estiver aplicado, nao faz nada.

Ref: PIPELINE_DOCS.md secao 10.3.
"""
import importlib.util
import os
import sys

# Shim global: relaxa o coerce de dtype do pandas (inserido apos 'import pandas as pd')
SHIM = '''# pandas>=2 compatibility shim (phages-pipeline): allow mixed-dtype .loc assignments
try:
    import pandas.core.internals.blocks as _pd_blk
    _orig_coerce = _pd_blk.Block.coerce_to_target_dtype
    def _patched_coerce(self, element, raise_on_upcast=False, *a, **k):
        return _orig_coerce(self, element, raise_on_upcast=False)
    _pd_blk.Block.coerce_to_target_dtype = _patched_coerce
except Exception:
    pass
'''

# Por arquivo: (anchor da 1a ocorrencia, linha-fix a inserir logo apos, assinatura idempotencia)
ASTYPE_PATCHES = {
    "phagcn.py": (
        "df = df.reset_index(drop=True)",
        'df["cluster"] = df["cluster"].astype(object)  # pandas>=2 fix (phages-pipeline)',
        'df["cluster"] = df["cluster"].astype(object)',
    ),
    "cherry.py": (
        "df = df.reset_index(drop=True)",
        'df["Score"] = df["Score"].astype(object)  # pandas>=2 fix (phages-pipeline)',
        'df["Score"] = df["Score"].astype(object)',
    ),
    "votu.py": (
        "cluster_df.reset_index(drop=True, inplace=True)",
        'cluster_df["cluster"] = cluster_df["cluster"].astype(object)  # pandas>=2 fix (phages-pipeline)',
        'cluster_df["cluster"] = cluster_df["cluster"].astype(object)',
    ),
}


def find_pkg_dir():
    spec = importlib.util.find_spec("phabox2")
    if spec and spec.submodule_search_locations:
        return list(spec.submodule_search_locations)[0]
    return None


def patch_shim(pkg_dir):
    path = os.path.join(pkg_dir, "phabox2.py")
    if not os.path.exists(path):
        return f"  phabox2.py: NAO encontrado (pulado)"
    with open(path, encoding="utf-8") as fh:
        src = fh.read()
    if "_patched_coerce" in src:
        return "  phabox2.py: shim ja aplicado"
    lines = src.splitlines(keepends=True)
    out, done = [], False
    for ln in lines:
        out.append(ln)
        if not done and ln.strip() == "import pandas as pd":
            out.append(SHIM)
            done = True
    if not done:
        return "  phabox2.py: anchor 'import pandas as pd' NAO encontrado (shim NAO aplicado)"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("".join(out))
    return "  phabox2.py: shim APLICADO"


def patch_astype(pkg_dir, fname, anchor, fix_line, signature):
    path = os.path.join(pkg_dir, fname)
    if not os.path.exists(path):
        return f"  {fname}: NAO encontrado (pulado)"
    with open(path, encoding="utf-8") as fh:
        src = fh.read()
    if signature in src:
        return f"  {fname}: ja aplicado"
    lines = src.splitlines(keepends=True)
    out, done = [], False
    for ln in lines:
        out.append(ln)
        if not done and anchor in ln:
            indent = ln[: len(ln) - len(ln.lstrip())]
            nl = "\n" if ln.endswith("\n") else ""
            out.append(f"{indent}{fix_line}{nl}")
            done = True
    if not done:
        return f"  {fname}: anchor NAO encontrado (NAO aplicado) -- verifique a versao do PhaBox2"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("".join(out))
    return f"  {fname}: APLICADO"


def main():
    pkg_dir = find_pkg_dir()
    if not pkg_dir:
        print("[patch_phabox2] phabox2 nao importavel neste env — nada a fazer.")
        return 0
    print(f"[patch_phabox2] pacote: {pkg_dir}")
    print(patch_shim(pkg_dir))
    for fname, (anchor, fix_line, sig) in ASTYPE_PATCHES.items():
        print(patch_astype(pkg_dir, fname, anchor, fix_line, sig))
    print("[patch_phabox2] OK.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
