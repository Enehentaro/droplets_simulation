# Droplets Simulation
Simulation of Virus-Laden Droplets Behavior in AFDET

## ドキュメント
https://enehentaro.github.io/droplets_simulation/

※[FORD](https://github.com/Fortran-FOSS-Programmers/ford)を使用しています.

## 使い方
  
### 依存関係解決・コンパイル&実行
  - ビルドに`cmake`コマンドと`make`コマンドを使うので, 要インストール.
  - [CMakeに関する説明はこちら](https://qiita.com/ijknabla/items/05270ae5e597705d0dae#cmake-%E3%81%AE%E3%82%A4%E3%83%B3%E3%82%B9%E3%83%88%E3%83%BC%E3%83%AB)

  ルートディレクトリ(README.mdのある階層)での作業.
  1. 「SampleCase」ディレクトリを複製したのち, 名前を変更する. (ケース名を付ける)
  1. ケースディレクトリ内の条件ファイル(condition.nml, initial_position.csv)を編集する.
  1. 以下OSに応じてコンパイル&実行する.

  **Linuxの場合**
  
  - `$ source build.sh`でビルド・CTest・パス通しまで出来る. (コンパイル手法を変える場合は各自で編集)
  - `$ MAIN`で実行. (プログラム名は自分で変える)

  **Windowsの場合**
  - `> .\build.bat`(Windows powershellにて実行)でビルドする. 
    - batファイル内の`cmake .. `の後に`-G "MinGW Makefiles"`でGeneratorを指定する.
  - `> .\build\bin\main.exe`で実行. (プログラム名は自分で変える)
  ※おそらくVSCodeでコンパイラ指定のポップアップウィンドウが出るので, そこでgfortranのパスを通す.
  
### シェルスクリプトの説明
  - シェルスクリプトを使用しない場合は, 下記順序に従ってコマンドを入力

  1. `$ mkdir build`でビルドディレクトリ作成する. 
  1. `$ cd build`で移動する. 
  1. `$ cmake ..`で依存関係解決.
      - `-D CMAKE_Fortran_COMPILER=[ifort/gfortran]`でコンパイラ指定.
      - `-D CMAKE_BUILD_TYPE=Debug`でデバッグ用コンパイルオプション付与.
      - `-D use_OpenMP=ON`でOpenMP用コンパイルオプション付与.
        - `find_package`でエラーが起こる場合はコンパイラにOpenMPが付属していないので, 要インストール. (tdm-gccはデフォルトでついてない)
  1. `$ make`でコンパイル.

## 条件ファイル(condition.nml, initial_position.csv)解説
### condition.nml
  - **リスタート位置 num_restart**
    - 通常は`0`を指定する.
    - `1以上`にすると, その値に対応するbackupファイルが読み込まれ, そこからリスタートが始まる.
  - **初期分布ファイル名 initialDistributionFName**
    - 指定したbackupファイル(.bu)が読み込まれ, それを飛沫初期分布とする.
    - 初期分布を固定したくない場合はコメントアウトする.
  - **飛沫周期発生 periodicGeneration**
    - 1秒当たりの発生飛沫数(整数)を指定する.
    - 初期配置飛沫をすべてNonActiveにしたのち, 順次Activateしていくので, 初期配置数が飛沫数の上限となる.
  - **気流データファイル名 path2FlowFile**
    - 実行ディレクトリからの相対パス, もしくは絶対パスを指定する.
    - 現在可能な流れ場ファイル:
      - VTK
      - INP
      - FLD
      - FPH
    - CUBE格子(PLOT3D)は、予め非構造格子に変換してから計算してください.
    - .arrayファイルを指定する場合、別途メッシュファイルが必要なので,`meshFile = ***`と指定する.
  - **ステップ数オフセット OFFSET**
    - 飛沫計算を, 流体連番ファイルの途中の番号から始めたいときに指定する.
  - **気流データを周期的に用いる場合の先頭と末尾 LoopHead, LoopTail**
    - 任意の区間の流体連番ファイルを繰り返し用いるときに指定する. (例えば呼吸のサイクル)
    - `(先頭) = (末尾)` とすると, そのステップ数到達後は流れ場の更新が起こらなくなる.
    - `(先頭) > (末尾)` とすれば, 特殊な処理は起こらず, 流体連番ファイルが順番に読み込まれる.
### initial_position.csv
  - 初期飛沫の配置帯（直方体）を設定する.
  - 左から順に, 直方体の中心座標(x,y,z), 直方体の幅(x,y,z).
  - 改行すれば配置帯を複数設定できる.


## 方程式

### 飛沫の蒸発方程式

  $$ \frac{dr}{dt} \space = \space -\left(1-\frac{RH}{100}\right) \cdot \frac{D e_{s}(T)}{\rho_{w} R_{v} T} \cdot \frac{1}{r} $$
  
  プログラム内では,2次精度ルンゲクッタ法で解いている.
  
### 飛沫の運動方程式

$$ m \frac{d \mathbf{v}}{dt} \space = \space m \mathbf{g} \space + \space C_D (\mathbf{v}) \space \cdot \space \frac{1}{2} \rho_a S \left | \mathbf{u}_a - \mathbf{v} \right | (\mathbf{u}_a - \mathbf{v}) $$

  プログラム内では,上式を無次元化・離散化した次式を解いている.
    
$$ \bar{\mathbf{v}}^{n + 1} \space = \space \frac{\bar{\mathbf{v}}^{n} \space + \space (\bar{\mathbf{g}} \space + \space C \bar{\mathbf{u}}_a)\Delta \bar{t}}{1 \space + \space C\Delta \bar{t}} \quad \left ( C \space = \space \frac{3 \rho_a}{8 \rho_w} \frac{C_D ( \mathbf{v}^{n} ) \left | \bar{\mathbf{u}}_a - \bar{\mathbf{v}}^{n} \right |}{\bar{r}^{n+1}} \right ) $$

## サブプログラム
  - CUBE2USG
    - CUBE格子を非構造格子に変換できる.
  - droplet2CSV
    - 飛沫計算結果を再度読み込み, 統計データ(浮遊数推移など)をCSVファイルに書き出す.
  - dropletCount
    - 飛沫計算結果を再度読み込み, カウントボックスを通過した飛沫数を調べる. optionディレクトリ内の"boxList.csv"を, ケースディレクトリに配置する必要がある.
  - initialTranslate
    - 飛沫の初期配置データを読み込み, 任意の座標への回転, 平行移動を行う. **by Konishi**

## おまけ機能
  - **複数ケース連続実行**
    - 実行時にTXTファイル名を入力すると, そのファイルに列挙された複数ケースを連続実行できる.
  - **basicSetting.nml**
    - optionディレクトリ内にある, 付着判定のオンオフや, 飛沫間合体の設定が可能.初期半径分布ファイルの指定も可能.
    
## CTest
  - コンパイル後, `$ ctest`でCTestが実行可. (buildディレクトリにて)
  - CTestの実行ディレクトリは、`test/`になる. (buildディレクトリではない)
  - テスト用プログラムはすべて`test/`で管理しよう.
