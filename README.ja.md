# 🌀 spintronics

純粋なRustで実装されたスピントロニクス・磁気ダイナミクス・トポロジカル物質シミュレーションライブラリ

**齊藤英治教授グループ (東京大学 / 理研CEMS) の先駆的研究成果に基づく実装**

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)]()
[![License](https://img.shields.io/badge/license-MIT%2FApache--2.0-blue)]()
[![Rust Version](https://img.shields.io/badge/rust-2021-orange)]()

## 🚀 概要

`spintronics` は、スピントロニクスおよび量子物性現象のシミュレーションを行うための包括的なRustクレートです。`scirs2` 科学計算エコシステム上に構築され、Rustの型安全性とゼロコスト抽象化を活用して、高速・安全・物理的に正確なシミュレーションを実現します：

- **スピンポンピング・輸送**: スピン流の生成と伝播
- **スピン・電荷変換**: 逆スピンホール効果 (ISHE)、スピンゼーベック効果 (SSE)
- **磁化ダイナミクス**: Landau-Lifshitz-Gilbert (LLG) 方程式ソルバー
- **トポロジカル現象**: スキルミオン、磁壁、トポロジカル電荷
- **ナノメカニカル結合**: Barnett効果、Einstein-de Haas効果
- **物理リザバー計算**: マグノンベースのニューロモルフィック計算
- **空洞マグノニクス**: マグノン-光子ハイブリッド系

*「Pythonのループは遅すぎる。でもC++のメモリ管理は疲れる」* - このライブラリは、コンパイル言語のパフォーマンスとRustの安全性を両立したい研究者・学生のために設計されています。

## 📊 開発状況

**現在のバージョン**: 0.2.0 ✅ **本番環境対応**

**最新リリース**: 2025年12月

### バージョン 0.2.0 ハイライト
- ✅ **インタラクティブWebデモ**: 4つの物理シミュレーションを備えたHTMX + Axumデモサブクレート
- ✅ **Pythonバインディング (PyO3)**: Pythonからネイティブパフォーマンスで利用可能
- ✅ **HDF5エクスポート**: 大規模データセット用のストレージ
- ✅ **メモリプールアロケータ**: ホットパスでのアロケーションを99%削減
- ✅ **Serdeシリアライゼーション**: JSON/バイナリデータ交換
- ✅ **単位検証**: 物理量の妥当性ランタイムチェック
- ✅ **パフォーマンス**: ホットパス関数に21個のinline属性
- ✅ **17のサンプル**: 難易度別に整理（初級/中級/上級）
- ✅ **448テスト合格**: 431ライブラリ (381 unit + 50 doc) + 17デモ、警告ゼロ

### コア機能
- ✅ **18実装モジュール**: 基礎物理から先端現象まで包括的カバー
- ✅ **60+ソースファイル**: よく整理されたモジュラー設計
- ✅ **5つの実験検証**: 著名論文との照合（Saitoh 2006, Woo 2016等）
- ✅ **インタラクティブWebデモ**: モダンなHTMX + Axumサブクレート
- ✅ **WebAssemblyサポート**: ブラウザベースシミュレーション対応
- ✅ **マルチプラットフォームCI/CD**: Ubuntu、macOS、Windowsでテスト済み
- ✅ **本番品質**: 警告ゼロ、448テスト合格

## ✨ 主要機能

### パフォーマンスと安全性
- ⚡ **高性能**: 純粋なRustで最適化された数値カーネル、SIMD対応
- 🛡️ **型安全**: Rustの所有権システムがスピン/角運動量の「消失」をコンパイル時に防止
- 🔒 **メモリ安全**: セグフォなし、データ競合なし、未定義動作なし
- 🎯 **ゼロコスト抽象化**: 物理的抽象が効率的な機械語にコンパイル

### 科学計算
- 📚 **物理学に整合したアーキテクチャ**: コード構造がハミルトニアンや輸送方程式に直接対応
- 🧮 **検証済みモデル**: 査読済み実験論文に基づく実装
- 📊 **再現可能な結果**: 制御された乱数シードによる決定論的シミュレーション
- 🔬 **実験検証**: 公表された実験結果を再現するサンプル

### 開発者体験
- 📖 **充実したドキュメント**: LaTeX数式を含む包括的なdocコメント
- 🧪 **徹底したテスト**: 物理的正しさを検証するユニットテスト
- 🔧 **最小限の依存**: 高速コンパイル、容易な統合
- 🌐 **エコシステム統合**: `scirs2` 科学計算スイートの一部
- 🐍 **Python統合**: PyO3バインディングでシームレスなPython連携
- 💾 **データエクスポート**: HDF5、JSON、CSV、VTK形式サポート
- ✅ **単位検証**: 物理量の妥当性チェック用14バリデーター

## 📚 主要参考文献

- E. Saitoh et al., "Conversion of spin current into charge current at room temperature: Inverse spin-Hall effect", *Appl. Phys. Lett.* **88**, 182509 (2006)
- K. Uchida et al., "Observation of the spin Seebeck effect", *Nature* **455**, 778-781 (2008)

## 📦 実装済みモジュール

18の物理学重視モジュール：

| モジュール | 物理概念 | 主要論文 / 概念 |
|--------|----------------|----------------------|
| **constants** | 物理定数 | ℏ, γ, e, μ_B, k_B、20+のNIST検証済み定数 |
| **vector3** | 3Dベクトル演算 | スピン/磁化演算に最適化 |
| **material** | 物質物性 | 強磁性体（YIG, Py）、界面、2D材料、トポロジカル材料 |
| **dynamics** | 磁化ダイナミクス | RK4、Heun、適応的手法を用いたLLGソルバー |
| **transport** | スピン輸送 | スピンポンピング（Saitoh 2006）、拡散方程式 |
| **effect** | スピン・電荷変換 | ISHE、SSE、SOT、Rashba、トポロジカルホール |
| **magnon** | マグノン伝播 | スピン波ダイナミクス、スピン鎖、マグノン検出 |
| **thermo** | 熱電効果 | 異常ネルンスト、熱マグノン、多層膜 |
| **texture** | 磁気テクスチャ | スキルミオン、磁壁、DMI、トポロジカル電荷 |
| **circuit** | スピン回路理論 | 抵抗ネットワーク、スピン蓄積 |
| **fluid** | スピン-渦結合 | 液体金属中のBarnett効果 |
| **mech** | ナノメカニカルスピントロニクス | Barnett、Einstein-de Haas、カンチレバー結合 |
| **ai** | 物理リザバー計算 | ニューロモルフィック計算用マグノンダイナミクス |
| **afm** | 反強磁性ダイナミクス | THzスピントロニクス（NiO, MnF₂等） |
| **stochastic** | 熱揺らぎ | 有限温度効果、Langevinダイナミクス |
| **cavity** | 空洞マグノニクス | マグノン-光子ハイブリッド量子系 |
| **memory** | メモリ管理 | プールアロケーター、ワークスペースバッファ (v0.2.0) |
| **units** | 単位検証 | 物理量の14バリデーター (v0.2.0) |
| **visualization** | データエクスポート | HDF5、JSON、CSV、VTK形式 (v0.2.0) |
| **python** | Pythonバインディング | Pythonユーザー向けPyO3統合 (v0.2.0) |

## クイックスタート

```rust
use spintronics::prelude::*;

// 物質のセットアップ (YIG/Pt系)
let yig = Ferromagnet::yig();
let interface = SpinInterface::yig_pt();
let pt_strip = InverseSpinHall::platinum();

// 磁化状態の初期化
let m = Vector3::new(1.0, 0.0, 0.0);
let h_ext = Vector3::new(0.0, 0.0, 1.0);

// LLG方程式を解く
let dm_dt = calc_dm_dt(m, h_ext, GAMMA, yig.alpha);

// スピンポンピング電流を計算
let js = spin_pumping_current(&interface, m, dm_dt);

// ISHEで電場に変換
let e_field = pt_strip.convert(interface.normal, js);
```

## 🎯 インストール

`Cargo.toml` に追加:

```toml
[dependencies]
spintronics = "0.2.0"
```

### オプション機能

```toml
[dependencies]
spintronics = { version = "0.2.0", features = ["python", "hdf5", "serde"] }
```

利用可能な機能:
- `python` - PyO3を使ったPythonバインディング
- `hdf5` - HDF5ファイルエクスポートサポート
- `serde` - JSON/バイナリシリアライゼーション
- `fem` - 有限要素法ソルバー
- `wasm` - WebAssemblyサポート

リポジトリから直接インストール:

```bash
git clone https://github.com/cool-japan/spintronics.git
cd spintronics
cargo build --release
```

## 💡 サンプル

**17の包括的サンプル**を難易度別に提供。完全ガイドは [`examples/README.md`](examples/README.md) を参照。

### 📚 クイックスタートサンプル（初級）

#### 1. YIG/Pt スピンポンピング + ISHE
```bash
cargo run --release --example yig_pt_pumping
```
齊藤らの著名な実験（2006）を再現：
- YIGにおける強磁性共鳴
- スピンポンピングによるスピン流生成
- Pt中の逆スピンホール効果による電圧検出

### 🔬 中級サンプル

- **スキルミオンダイナミクス** - トポロジカルスピンテクスチャと電流駆動運動
- **スピントルクオシレーター** - 自励振動と位相同期
- **2D材料** - ファンデルワールスヘテロ構造（CrI₃, Fe₃GeTe₂）

### 🚀 上級サンプル

- **FEMマイクロマグネティクス** - 現実的な形状の有限要素シミュレーション
- **並列マグノンダイナミクス** - マルチスレッド大規模シミュレーション
- **熱マグノン輸送** - スピンゼーベック効果と熱勾配
- **トポロジカル絶縁体** - 表面状態とEdelstein効果
- **リザバー計算** - マグノンを使ったニューロモルフィック計算

詳細は [`examples/README.md`](examples/README.md) を参照：
- 全17サンプルの詳細説明
- 背景知識別の学習パス
- 難易度評価と前提知識
- 機能要件とビルドコマンド
- リザバー計算の物理実装

### 全サンプルの実行
```bash
# サンプルを順次実行
for example in yig_pt_pumping magnon_propagation advanced_spintronics \
               mech_coupling fluid_barnett reservoir_computing; do
    cargo run --release --example $example
done
```

## 🌐 インタラクティブWebデモ

**v0.2.0の新機能！** インタラクティブなWebデモを試す：

```bash
cd demo
cargo run --release
# ブラウザで http://localhost:3000 を開く
```

### 利用可能なデモ

1. **LLG磁化ダイナミクス** (`/llg`)
   - リアルタイムLLGソルバーと軌道可視化
   - インタラクティブなパラメータ制御（減衰、磁場、初期状態）
   - 時間ステップ設定可能なRK4積分

2. **スピンポンピング計算機** (`/spin-pumping`)
   - Saitoh 2006 APL実験を再現
   - 物質選択（YIG、パーマロイ、CoFeB）
   - 周波数とRF磁場のパラメータスイープ

3. **材料エクスプローラー** (`/materials`)
   - 強磁性体間の磁気物性比較
   - 飽和磁化、減衰、交換スティフネス
   - 一般的なスピントロニクス材料データベース

4. **スキルミオン可視化** (`/skyrmion`)
   - リアルタイム磁化場レンダリング
   - ヘリシティ（Néel/Bloch）とカイラリティ制御
   - トポロジカル電荷計算

**技術スタック**: Axum + HTMX + Askama（サーバーサイドレンダリング、JavaScriptフレームワークなし）

詳細は [`demo/README.md`](demo/README.md) と [`demo/TESTING.md`](demo/TESTING.md) を参照。

## 🧪 テスト

全テストスイートを実行：

```bash
cargo test --all
```

出力付きでテスト実行：
```bash
cargo test -- --nocapture
```

特定のモジュールのテスト：
```bash
cargo test dynamics::
cargo test transport::
```

デモサブクレートのテスト：
```bash
cd demo
./test_server.sh  # 自動エンドポイントテスト
cargo test        # ユニットテスト
```

### テストカバレッジ

**合計: 448テスト合格** (431ライブラリ + 17デモ)
- ✅ **381ユニットテスト**: コア物理計算
- ✅ **50ドキュメントテスト**: ドキュメント内サンプル
- ✅ **17デモテスト**: Webエンドポイントと物理検証

全モジュールに包括的なテストを実装：
- ✅ **物理的正しさ**: 保存則、対称性、ゲージ不変性
- ✅ **エッジケース**: ゼロ磁場、平行/反平行配置、境界条件
- ✅ **物質パラメータ**: 文献値との検証
- ✅ **数値安定性**: 収束テスト、安定性解析
- ✅ **統合テスト**: 複数モジュール間の物理ワークフロー

## ⚡ パフォーマンス

Rustのゼロコスト抽象化とコンパイル時最適化により、インタープリタ言語に対して大幅なパフォーマンス向上を実現：

| ベンチマーク | Python (NumPy) | Rust (このライブラリ) | スピードアップ |
|----------|----------------|----------------------|---------------|
| LLG 1ステップ (1000スピン) | ~50 ms | ~0.8 ms | **62x** |
| スピンポンピング計算 | ~10 ms | ~0.05 ms | **200x** |
| マグノン伝播 (10,000ステップ) | ~5 s | ~80 ms | **62x** |

*ベンチマーク環境: Intel Core i7-10700K @ 3.8 GHz, 32 GB RAM*

### 最適化技術
- **SIMD自動ベクトル化**: コンパイラによるベクトル演算最適化
- **インライン展開**: ホットパス関数に21個の`#[inline]`属性
- **メモリプールアロケーター**: ホットパスでのアロケーションを99%削減
- **ゼロコストラッパー**: newtype パターンによる型安全性とランタイムオーバーヘッドゼロ

## 🤝 コントリビューション

貢献を歓迎します！ [`CONTRIBUTING.md`](CONTRIBUTING.md) を参照してください。

### 行動規範

本プロジェクトは [Contributor Covenant](CODE_OF_CONDUCT.md) 行動規範を採用しています。

### 開発ワークフロー

1. リポジトリをフォーク
2. 機能ブランチを作成 (`git checkout -b feature/amazing-physics`)
3. 変更をコミット (`git commit -m 'Add skyrmion Hall effect'`)
4. ブランチにプッシュ (`git push origin feature/amazing-physics`)
5. プルリクエストを開く

## 📄 ライセンス

このプロジェクトは MIT OR Apache-2.0 のデュアルライセンスです。

## 📧 コンタクト

- **メンテナー**: COOLJAPAN OÜ (Team KitaSan)
- **リポジトリ**: https://github.com/cool-japan/spintronics
- **ドキュメント**: https://docs.rs/spintronics
- **Issues**: https://github.com/cool-japan/spintronics/issues

## 🙏 謝辞

本ライブラリは以下の研究グループと論文に基づいています：

- **齊藤英治グループ** (東京大学 / 理研CEMS) - スピンポンピング、ISHE、SSEの先駆的研究
- **内田健一, 齊藤英治 et al.** - スピンゼーベック効果の発見 (Nature 2008)
- **Woo et al.** - 室温スキルミオンの観測 (Nature Materials 2016)

そして、世界中のスピントロニクス研究コミュニティに感謝します。

---

**Copyright © 2025 COOLJAPAN OÜ (Team KitaSan)**

Licensed under MIT OR Apache-2.0
