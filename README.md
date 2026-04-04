# Quantum-chemical-calculations

Pythonを用いた量子化学計算（制限ハートリー・フォック法: RHF）および分子の構造最適化を行うプログラム群です。量子化学計算の基礎的なアルゴリズムの実装・学習・検証を目的としています。

## 特徴

- **RHF法によるエネルギー計算**: Roothaan方程式を反復的に解き、分子の全エネルギーを計算します。
- **構造最適化 (1次元 / 2次元)**: HeH+のような二原子分子の1次元構造最適化や、H2Oなどの2次元構造最適化をスクリプトから簡単に実行できます。
- **結果の可視化とデータ出力**: 最適化過程のエネルギー変化や座標を `pandas` を用いてCSVファイルに保存するほか、エネルギー曲線や2次元のエネルギー曲面のグラフ化をサポートしています。
- **汎用的なデータ形式への対応**: 分子の初期座標はXYZファイル (`.xyz`) から取得し、基底関数のデータ（STO-3Gなど）はJSONファイルから柔軟に読み込めるよう設計しています。

## 主要なファイルとディレクトリ

- **`RHF_s_orbital_only.py`**: s軌道のみを考慮したRHF法の計算スクリプトです。重なり積分、運動エネルギー、核引力、電子間反発積分（ERI）からFock行列を構築し、自己無撞着場（SCF）計算の具体的な実装を確認できます。
- **`1d_optimize_HeH+.py`**: HeH+の1次元構造最適化を実行し、エネルギー曲線のグラフを出力するスクリプトです。
- **`2d_optimize_H2O.py`**: H2Oの2次元構造最適化を実行し、実行時間の計測、CSVへの結果(`Structural_optimization of H2O.csv`)の保存、およびエネルギー曲面のグラフ化を行うスクリプトです。
- **`test.py`**: 各種モジュール（`libs.RHF`, `libs.optimize`）の動作確認用のテストコード群です。
- **`libs/`**: 計算のコアとなる自作モジュールが格納されています。
  - **`libs.RHF.RHF()`**: 分子のXYZ座標ファイルと基底関数のJSONファイルを読み込み、自己無撞着場（SCF）計算によって全エネルギーを求めるRHF計算のメイン関数です。
    - **引数**: `input_file` (str: XYZファイルのパス), `base_function_file` (str: 基底関数のJSONファイルのパス), `eps` (float: 収束閾値), `max_iter` (int: Roothaan方程式の繰り返し計算の最大回数), `print_E` (bool: 途中経過の表示フラグ) など。
    - **返り値**: 収束した分子の全エネルギー (float)。
  - **`libs.optimize.optimize()`**: 指定した変数の範囲とステップ幅で1次元のグリッドサーチを行い、エネルギーが最小となる構造を探索します（例: 二原子分子の結合長）。
    - **引数**: `raw_data` (str: 初期座標のXYZファイルパス), `steps` (list[float]: 探索のステップ幅), `val_ranges` (list[list[float]]: 探索する変数の範囲) とSCF計算用の各種パラメータ。
    - **返り値**: 探索した各点の座標変数とエネルギーを格納した辞書のリスト (`list[dict]`)。
  - **`libs.optimize.optimize_H2O_type()`**: 2つの変数を同時に変化させて2次元のグリッドサーチを行い、エネルギー最小構造を探索します（例: H2Oの結合長と結合角の最適化）。
    - **引数**: `optimize()` と同様ですが、`steps` や `val_ranges` に2次元分のリストを指定します。
    - **返り値**: 探索結果のリスト (`list[dict]`)、最小エネルギー値 (`min_E` : float)、エネルギーを最小とする変数 (`min_x` : float, `min_y` : float)。
  - **`libs.integrals` モジュール**: `nuclear_attraction`（核引力積分）や `electron_repulsion`（電子間反発積分）など、Boys関数を用いた各種分子積分を計算する関数が含まれています。
    - ※ この実装部分は、jjgoings/McMurchie-Davidson から引用しています。
    - ガウス軌道の軌道指数、角運動量、中心座標の各種パラメータを引数に受け取り、解析的な積分値 (float) を返します。
  - **`libs.optimize` のグラフ化関数**: `single_d_figure_make()` や `double_d_figure_make()` を用い、最適化過程のエネルギー曲線/曲面の可視化を行います。
    - `optimize()` などの返り値である結果のリストや、最小エネルギー時の座標を引数に取り、`matplotlib` によるグラフを描画します。
- **`samples/`**: 計算対象となる分子のXYZファイル（`H2O.xyz`, `HeH+_opt.xyz` など）が格納されています。
- **`base_function_data/`**: ガウス型基底関数（GTO）の係数と指数をまとめたJSONファイル（`STO-3G.json` など）が格納されています。

## 必須要件（Requirements）

実行には以下のPythonライブラリが必要です。
- `numpy`
- `scipy`
- `pandas`
- `matplotlib` （グラフ描画用）

```bash
pip install numpy scipy pandas matplotlib
```

## 使い方（Usage）

各スクリプトを直接実行することで、計算や最適化処理が走ります。

**例1: HeH+の1次元構造最適化を実行し、エネルギー曲線を出力する。**
```bash
python 1d_optimize_HeH+.py
```

**例2: H2Oの2次元構造最適化を実行し、エネルギー曲面のグラフを出力する。また、計算結果を格納したCSVファイルを出力する**
```bash
python 2d_optimize_H2O.py
```

**例3: s軌道のみを考慮したRHF法（基底関数の展開など）をテストする**
```bash
python RHF_s_orbital_only.py
```