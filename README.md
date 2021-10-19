# Droplets Simulation
Simulation of Droplets Behavior in AFDET

## main.f90
  飛沫計算メインプログラム。前処理プログラムは取り込まれました。
  実行時にファイルの有無から判断して、必要であれば前処理が勝手に行われます。
  現在可能な流れ場ファイル：
  - VTK
  - INP
  - FLD
  - PLOT3D

## 使い方
  ※環境はIntel Fortran、Linuxを想定しています。その他の環境では適宜書き換えが必要です。
  コンパイルに`make`コマンドを使います（makeのインストールが必要）。
  - Makefileのあるディレクトリで `make` コマンド（コンパイル）
  - 「sample_case」ディレクトリを複製したのち、名前を変更する（ケース名を付ける）
  - ケースディレクトリ内の条件ファイルを編集
  - 「droplet」という実行ファイルを実行し、ケース名を入力して計算開始

## 外部サブルーチン「management_droplet」
  dropletManager.f90内で定義されているサブルーチン「management_droplet」は、毎ステップ呼び出される外部サブルーチンです。
  自由に処理を追加することができます（例えば任意の範囲内にいる飛沫数のカウントなど）。 ご利用ください。

## 方程式

  解くべき方程式は次の通り。  
<img src="https://latex.codecogs.com/gif.latex?m&space;\frac{d&space;\mathbf{v}}{dt}&space;=&space;m&space;\mathbf{g}&space;&plus;&space;C_D&space;\cdot&space;\frac{1}{2}\rho_a&space;S&space;\left&space;|&space;\mathbf{u}_a&space;-&space;\mathbf{v}&space;\right&space;|(\mathbf{u}_a&space;-&space;\mathbf{v})" />

  プログラム内では、上式を無次元化・離散化した次式を解いている。  
<img src="https://latex.codecogs.com/gif.latex?\bar{\mathbf{v}}^{n&plus;1}&space;=&space;\frac{\bar{\mathbf{v}}^{n}&space;&plus;&space;(\bar{\mathbf{g}}&space;&plus;&space;C\bar{\mathbf{u}}_a)\Delta&space;\bar{t}}{1&plus;C\Delta&space;\bar{t}}" />
