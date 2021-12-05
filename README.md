# Droplets Simulation
Simulation of Virus Droplets Behavior in AFDET

## 使い方
  ※環境は **Intel Fortran, Linux** を想定しています。その他の環境では適宜書き換えが必要です。
  コンパイルに`make`コマンドを使います（makeのインストールが必要）。
  1. 「SampleCase」ディレクトリを複製したのち、名前を変更する（ケース名を付ける）。
  2. ケースディレクトリ内の条件ファイル(condition.txt, initial_position.csv)を編集。
  3. Makefileのあるディレクトリで `make` コマンド（コンパイル）。
  4. `./droplet`で実行。ケース名を入力して計算開始。

## 条件ファイル(condition.txt)解説
  - **リスタート位置**
    - 通常は`0`を指定
    - `1以上`にすると、その値に対応するbackupファイルが読み込まれ、そこからリスタートが始まる
    - `-1`にすると、InitialDistribution.buという名前のファイルが読み込まれ、それを飛沫初期分布とする。全く同じ初期分布から計算を始めたいときに使う。なお、InitialDistribution.buは、通常実行時にbackupディレクトリに最初に出力されるので、それをケースディレクトリに配置する必要がある。
  - **気流データファイル名**
    - 実行ディレクトリからのパスを指定（ケースディレクトリ基準ではないので注意）（改善の余地あり）
    - 現在可能な流れ場ファイル：
      - VTK
      - INP
      - FLD
    - CUBE格子(PLOT3D)は、予め非構造格子に変換してから計算してください。
  - **ステップ数オフセット**
    - 飛沫計算を、流体連番ファイルの途中の番号から始めたいときに指定
  - **気流データを周期的に用いる場合の先頭と末尾**
    - 任意の区間の流体連番ファイルを繰り返し用いるときに指定（例えば呼吸のサイクル）
    - `(末尾)＝(先頭)` とすると、そのステップ数到達後は流れ場の更新が起こらなくなる
    - `(末尾)＜(先頭)` とすれば、特殊な処理は起こらず、流体連番ファイルが順番に読み込まれる

## main.f90
  飛沫計算メインプログラム。前処理プログラムは取り込まれました。実行時にファイルの有無から判断して、必要であれば前処理が勝手に行われます。

## 外部サブルーチン「dropletManagement」
  dropletManager.f90内で定義されているサブルーチン「dropletManagement」は、毎ステップ呼び出される外部サブルーチンです。
  自由に処理を追加することができます（例えば任意の範囲内にいる飛沫数のカウントなど）。 ご利用ください。

## 方程式

  解くべき方程式は次の通り。  
<img src="https://latex.codecogs.com/gif.latex?m&space;\frac{d&space;\mathbf{v}}{dt}&space;=&space;m&space;\mathbf{g}&space;&plus;&space;C_D&space;\cdot&space;\frac{1}{2}\rho_a&space;S&space;\left&space;|&space;\mathbf{u}_a&space;-&space;\mathbf{v}&space;\right&space;|(\mathbf{u}_a&space;-&space;\mathbf{v})" />

  プログラム内では、上式を無次元化・離散化した次式を解いている。  
<img src="https://latex.codecogs.com/gif.latex?\bar{\mathbf{v}}^{n&plus;1}&space;=&space;\frac{\bar{\mathbf{v}}^{n}&space;&plus;&space;(\bar{\mathbf{g}}&space;&plus;&space;C\bar{\mathbf{u}}_a)\Delta&space;\bar{t}}{1&plus;C\Delta&space;\bar{t}}" />

## おまけ機能(optionディレクトリ内のTXTファイル)
  - **case_list.txt**
    - ここにケース名を複数列挙し、実行時に`option/case_list.txt`と入力すると、複数ケースを連続実行できる
