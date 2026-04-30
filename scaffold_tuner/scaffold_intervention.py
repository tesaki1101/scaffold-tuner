# -*- coding: utf-8 -*-
import random
from rdkit import Chem
from rdkit.Chem import Lipinski, Descriptors, rdRGroupDecomposition
from rdkit.Chem.Scaffolds import MurckoScaffold


# ──────────────────────────────────────────────
# 1. 親分子 → Murcko コア（Rラベル付き）
# ──────────────────────────────────────────────
def make_murcko_core_with_rlabels(parent_smiles, max_cuts=None):
    """
    親分子SMILESから Murcko scaffold を作り、
    サイドチェーンとの境界に [*:1], [*:2], ... を付けたコアを返す。

    Returns: (parent_mol, scaffold_mol, core_with_labels)
    """
    parent_mol = Chem.MolFromSmiles(parent_smiles)
    if parent_mol is None:
        raise ValueError("無効な親分子SMILESです")

    scaffold_mol = MurckoScaffold.GetScaffoldForMol(parent_mol)
    if scaffold_mol is None:
        raise ValueError("Murcko scaffoldの生成に失敗しました")

    matches = parent_mol.GetSubstructMatches(scaffold_mol)
    if not matches:
        raise ValueError("Murcko scaffoldが親分子にマッチしませんでした")

    # scaffold に含まれる親分子の原子インデックスと、その逆引きマップ
    match = matches[0]
    core_atom_idxs = set(match)
    parent_to_core = {p: c for c, p in enumerate(match)}

    # コアとサイドチェーンをまたぐ結合（境界結合）を収集
    boundary_bonds = [
        b for b in parent_mol.GetBonds()
        if (b.GetBeginAtomIdx() in core_atom_idxs) ^ (b.GetEndAtomIdx() in core_atom_idxs)
    ]
    if max_cuts is not None:
        boundary_bonds = boundary_bonds[:max_cuts]

    # 境界結合ごとに、コア側原子へダミー原子 [*:N] を追加
    rw_core = Chem.RWMol(scaffold_mol)
    for rlabel, bond in enumerate(boundary_bonds, start=1):
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        core_parent_idx = a1 if a1 in core_atom_idxs else a2
        core_idx = parent_to_core[core_parent_idx]

        dummy = Chem.Atom(0)          # ダミー原子 *
        dummy.SetAtomMapNum(rlabel)
        dummy_idx = rw_core.AddAtom(dummy)
        rw_core.AddBond(core_idx, dummy_idx, bond.GetBondType())

    core_with_labels = rw_core.GetMol()
    Chem.SanitizeMol(core_with_labels)
    return parent_mol, scaffold_mol, core_with_labels


# ──────────────────────────────────────────────
# 2. ユーティリティ
# ──────────────────────────────────────────────
def count_features(mol):
    """HBA, HBD, 芳香環数, 回転可能結合数, ヘテロ原子数 をタプルで返す"""
    return (
        Lipinski.NumHAcceptors(mol),
        Lipinski.NumHDonors(mol),
        Descriptors.NumAromaticRings(mol),
        Descriptors.NumRotatableBonds(mol),
        Descriptors.NumHeteroatoms(mol),
    )


def get_r_labels(core_mol):
    """コア中のダミー原子 [*:N] の map番号リストを返す"""
    return sorted({
        atom.GetAtomMapNum()
        for atom in core_mol.GetAtoms()
        if atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() > 0
    })


# ──────────────────────────────────────────────
# 3. 親分子から元のRグループを取得
# ──────────────────────────────────────────────
def get_original_rgroups(parent_smiles, core_with_labels):
    """
    Rグループ分解により、親分子の各置換基を {"R1": SMILES, ...} で返す
    """
    parent_mol = Chem.MolFromSmiles(parent_smiles)
    if parent_mol is None:
        raise ValueError("無効な親分子SMILESです")

    params = rdRGroupDecomposition.RGroupDecompositionParameters()
    params.labels = rdRGroupDecomposition.RGroupLabels.AtomMapLabels
    params.rgroupLabelling = rdRGroupDecomposition.RGroupLabelling.AtomMap
    params.onlyMatchAtRGroups = True

    rgd = rdRGroupDecomposition.RGroupDecomposition(core_with_labels, params)
    if rgd.Add(parent_mol) < 0:
        raise ValueError("RGroupDecomposition.Add() が失敗しました")
    if not rgd.Process():
        raise RuntimeError("RGroupDecomposition.Process() が失敗しました")

    rows = rgd.GetRGroupsAsRows()
    if not rows:
        raise RuntimeError("Rグループ分解の結果が得られませんでした")

    # "Core" キーを除き、SMILESに変換して返す
    return {
        key: Chem.MolToSmiles(mol)
        for key, mol in rows[0].items()
        if key != "Core" and mol is not None
    }


# ──────────────────────────────────────────────
# 4. フラグメントの結合
# ──────────────────────────────────────────────
def _get_dummy_info(mol, map_num):
    """[*:map_num] のダミー原子を探し (dummy_idx, neighbor_idx, bond_type) を返す"""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() == map_num:
            nbr = atom.GetNeighbors()[0]
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            return atom.GetIdx(), nbr.GetIdx(), bond.GetBondType()
    raise ValueError(f"ダミー原子 [*:{map_num}] が見つかりません")


def attach_fragment(scaffold_mol, fragment_mol, map_num):
    """scaffold と fragment の [*:map_num] 同士を結合し、ダミー原子を除去する"""
    s_dummy, s_nbr, s_bond = _get_dummy_info(scaffold_mol, map_num)
    f_dummy, f_nbr, _      = _get_dummy_info(fragment_mol, map_num)

    offset = scaffold_mol.GetNumAtoms()
    emol = Chem.EditableMol(Chem.CombineMols(scaffold_mol, fragment_mol))
    emol.AddBond(s_nbr, f_nbr + offset, order=s_bond)

    # インデックスがずれないよう大きい方から削除
    for idx in sorted([s_dummy, f_dummy + offset], reverse=True):
        emol.RemoveAtom(idx)

    mol = emol.GetMol()
    Chem.SanitizeMol(mol)
    return mol


def build_molecule(scaffold_smiles, rgroups):
    """
    scaffold_smiles に対して rgroups（例: {"R1": "[*:1]N", "R2": "[*:2]OC"}）を順に結合する
    """
    mol = Chem.MolFromSmiles(scaffold_smiles)
    if mol is None:
        raise ValueError("無効なscaffold SMILESです")

    # R番号の昇順に結合（順序依存を避けるため）
    for site in sorted(rgroups, key=lambda x: int(x[1:])):
        frag = Chem.MolFromSmiles(rgroups[site])
        if frag is None:
            raise ValueError(f"無効なフラグメントSMILES: {rgroups[site]}")
        mol = attach_fragment(mol, frag, int(site[1:]))
    return mol


# ──────────────────────────────────────────────
# 5. 各特徴量を増減させる置換基プール
# ──────────────────────────────────────────────

# [*:1] をテンプレートとして定義し、使用時に map_num へ置き換える
_FRAGMENT_TEMPLATES = {
    # HBD（水素結合ドナー）
    "add_hbd": [
        "[*:1]N",           # アミノ基
        "[*:1]NC",          # メチルアミノ基
        "[*:1]O",           # ヒドロキシ基
        "[*:1]NS(=O)(=O)C", # スルホンアミド様
    ],
    "remove_hbd": [
        "[*:1]C",           # メチル基（HBDなし）
        "[*:1]CC",          # エチル基（HBDなし）
        "[*:1]F",           # フルオロ基（HBDなし）
        "[*:1]Cl",          # クロロ基（HBDなし）
        "[*:1]C(C)C",       # イソプロピル基（HBDなし）
    ],
    # HBA（水素結合アクセプター）
    "add_hba": [
        "[*:1]F",           # フルオロ基
        "[*:1]OC",          # メトキシ基
        "[*:1]C#N",         # シアノ基
        "[*:1]C(=O)C",      # ケトン様
        "[*:1]N(C)C",       # 三級アミン
    ],
    "remove_hba": [
        "[*:1]C",           # メチル基（HBAなし）
        "[*:1]CC",          # エチル基（HBAなし）
        "[*:1]F",           # フルオロ基（HBAとして数えられない場合が多い）
        "[*:1]Cl",          # クロロ基（HBAなし）
        "[*:1]c1ccccc1",    # フェニル基（HBAなし）
    ],
    # 芳香環数（NumAromaticRings）
    "add_ar": [
        "[*:1]c1ccccc1",    # フェニル基
        "[*:1]c1ccncc1",    # ピリジル基
        "[*:1]c1ccoc1",     # フラニル基
        "[*:1]c1ccsc1",     # チエニル基
    ],
    "remove_ar": [
        "[*:1]C",           # メチル基（芳香環なし）
        "[*:1]CC",          # エチル基（芳香環なし）
        "[*:1]CCC",         # プロピル基（芳香環なし）
        "[*:1]C(C)C",       # イソプロピル基（芳香環なし）
        "[*:1]C1CCCCC1",    # シクロヘキシル基（芳香環なし）
    ],
    # 回転可能結合数（NumRotatableBonds）
    "add_rb": [
        "[*:1]CC",          # エチル基（回転可能結合+1）
        "[*:1]CCC",         # プロピル基（回転可能結合+2）
        "[*:1]CCCC",        # ブチル基（回転可能結合+3）
        "[*:1]COC",         # メトキシメチル基（回転可能結合+2）
        "[*:1]CC(=O)C",     # アセトニル基（回転可能結合+2）
    ],
    "remove_rb": [
        "[*:1]C",           # メチル基（回転可能結合なし）
        "[*:1]F",           # フルオロ基（回転可能結合なし）
        "[*:1]Cl",          # クロロ基（回転可能結合なし）
        "[*:1]C1CC1",       # シクロプロピル基（環で回転制限）
        "[*:1]C1CCCCC1",    # シクロヘキシル基（環で回転制限）
    ],
    # ヘテロ原子数（NumHeteroatoms）
    "add_ha": [
        "[*:1]N",           # アミノ基（N追加）
        "[*:1]O",           # ヒドロキシ基（O追加）
        "[*:1]OC",          # メトキシ基（O追加）
        "[*:1]F",           # フルオロ基（F追加）
        "[*:1]S",           # チオール基（S追加）
    ],
    "remove_ha": [
        "[*:1]C",           # メチル基（ヘテロ原子なし）
        "[*:1]CC",          # エチル基（ヘテロ原子なし）
        "[*:1]C(C)C",       # イソプロピル基（ヘテロ原子なし）
        "[*:1]c1ccccc1",    # フェニル基（ヘテロ原子なし）
        "[*:1]C1CCCCC1",    # シクロヘキシル基（ヘテロ原子なし）
    ],
}

# mode → (特徴量インデックス, 増減方向)
# 特徴量インデックス: count_features の返り値に対応
#   0=HBA, 1=HBD, 2=NumAromaticRings, 3=NumRotatableBonds, 4=NumHeteroatoms
# 増減方向: +1=増加を期待, -1=減少を期待
_MODE_CONFIG = {
    "add_hbd":    (1, +1),
    "remove_hbd": (1, -1),
    "add_hba":    (0, +1),
    "remove_hba": (0, -1),
    "add_ar":     (2, +1),
    "remove_ar":  (2, -1),
    "add_rb":     (3, +1),
    "remove_rb":  (3, -1),
    "add_ha":     (4, +1),
    "remove_ha":  (4, -1),
}

def get_fragment_pool(mode, map_num):
    """mode に対応する置換基テンプレートを map_num 付きで返す"""
    if mode not in _FRAGMENT_TEMPLATES:
        raise ValueError(f"mode は {list(_MODE_CONFIG.keys())} のいずれかを指定してください")
    return [t.replace("[*:1]", f"[*:{map_num}]") for t in _FRAGMENT_TEMPLATES[mode]]


# ──────────────────────────────────────────────
# 6. メイン：候補構造の提案
# ──────────────────────────────────────────────
def propose_structures(
    parent_smiles,
    mode="add_hbd",
    max_cuts=None,
    max_candidates=10,
    shuffle_sites=False,
    shuffle_fragments=False,
    random_seed=None,
):
    """
    親分子の特徴量を1つ増減させた候補構造を提案する。
    変更しない部位には親分子の元の置換基をそのまま戻す。

    Parameters
    ----------
    parent_smiles    : 親分子のSMILES
    mode             : 以下のいずれかを指定する
                         "add_hbd"    : 水素結合ドナー数を増やす
                         "remove_hbd" : 水素結合ドナー数を減らす
                         "add_hba"    : 水素結合アクセプター数を増やす
                         "remove_hba" : 水素結合アクセプター数を減らす
                         "add_ar"     : 芳香環数を増やす
                         "remove_ar"  : 芳香環数を減らす
                         "add_rb"     : 回転可能結合数を増やす
                         "remove_rb"  : 回転可能結合数を減らす
                         "add_ha"     : ヘテロ原子数を増やす
                         "remove_ha"  : ヘテロ原子数を減らす
    max_cuts         : scaffoldに付与する置換可能部位数の上限（None=全部位）
    max_candidates   : 返す最大候補数
    shuffle_sites    : True のとき部位をランダム順で探索
    shuffle_fragments: True のとき置換基候補をランダム順で探索
    random_seed      : 乱数シード
    """
    if mode not in _MODE_CONFIG:
        raise ValueError(f"mode は {list(_MODE_CONFIG.keys())} のいずれかを指定してください")

    if random_seed is not None:
        random.seed(random_seed)

    parent_mol, _, core = make_murcko_core_with_rlabels(parent_smiles, max_cuts=max_cuts)
    core_smiles = Chem.MolToSmiles(core)
    labels = get_r_labels(core)
    if not labels:
        return []

    parent_canon = Chem.MolToSmiles(parent_mol)
    parent_feats = count_features(parent_mol)
    feat_idx, direction = _MODE_CONFIG[mode]  # 比較する特徴量と増減方向

    original_rgroups = get_original_rgroups(parent_smiles, core)

    if shuffle_sites:
        labels = random.sample(labels, len(labels))

    results, seen = [], set()

    for label in labels:
        site = f"R{label}"
        pool = get_fragment_pool(mode, label)
        if shuffle_fragments:
            pool = random.sample(pool, len(pool))

        for frag_smiles in pool:
            try:
                # 1箇所だけ置換基を差し替えて分子を構築
                rgroups = {**original_rgroups, site: frag_smiles}
                gen_mol = build_molecule(core_smiles, rgroups)
                gen_canon = Chem.MolToSmiles(gen_mol)

                # 親分子と同じ or 既出の場合はスキップ
                if gen_canon == parent_canon or gen_canon in seen:
                    continue
                seen.add(gen_canon)

                gen_feats = count_features(gen_mol)
                # 対象特徴量が direction の方向に変化しているかチェック
                if direction * gen_feats[feat_idx] > direction * parent_feats[feat_idx]:
                    results.append({
                        "mode": mode,
                        "site": site,
                        "fragment": frag_smiles,
                        "parent_smiles": parent_canon,
                        "generated_smiles": gen_canon,
                        "parent_features": parent_feats,
                        "generated_features": gen_feats,
                    })
            except Exception:
                continue

            if len(results) >= max_candidates:
                return results

    return results


# ──────────────────────────────────────────────
# 7. 結果の表示
# ──────────────────────────────────────────────
FEATURE_NAMES = ["NumHAcceptors", "NumHDonors", "NumAromaticRings", "NumRotatableBonds", "NumHeteroatoms"]
# count_features の返り値の順番と対応

def print_proposals(results):
    if not results:
        print("有効な候補が見つかりませんでした")
        return

    for i, rec in enumerate(results, start=1):
        print(f"=== Candidate {i} ===")
        print(f"mode: {rec['mode']},  site: {rec['site']},  fragment: {rec['fragment']}")
        print(f"parent SMILES   : {rec['parent_smiles']}")
        print(f"generated SMILES: {rec['generated_smiles']}")
        print(f"features {FEATURE_NAMES}")
        print(f"  parent   : {rec['parent_features']}")
        print(f"  generated: {rec['generated_features']}")
        print()
