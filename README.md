# Droplets Simulation
Simulation of Virus-Laden Droplets Behavior in AFDET

## ドキュメント
https://enehentaro.github.io/droplets_simulation/

※[FORD](https://github.com/Fortran-FOSS-Programmers/ford)を使用しています

## 使い方
  
### 依存関係解決・コンパイル
  - ビルドに`cmake`コマンドを使います（[CMakeのインストール](https://qiita.com/ijknabla/items/05270ae5e597705d0dae#cmake-%E3%81%AE%E3%82%A4%E3%83%B3%E3%82%B9%E3%83%88%E3%83%BC%E3%83%AB)が必要）。
  
  - `$ source build.sh`でビルド・CTest・パス通しまで出来る(コンパイル手法を変える場合は各自で編集)
  - `$ MAIN`で実行
  
#### `build.sh`を使わない場合（Windowsだとそもそもシェルスクリプト使えないかも）
  1. `$ mkdir build`でビルドディレクトリ作成
  1. `$ cd build`で移動
  1. `$ cmake ..`で依存関係解決
      - 一部の環境(Windows等)では、`-G "MinGW Makefiles"`でGenerator指定
      - `-D CMAKE_Fortran_COMPILER=[ifort/gfortran]`でコンパイラ指定
      - `-D CMAKE_BUILD_TYPE=debug`でデバッグ用コンパイルオプション付与
  1. `$ make`でコンパイル
  
#### 実行
  ルートディレクトリ（README.mdのあるディレクトリ）での作業
  
  1. 「SampleCase」ディレクトリを複製したのち、名前を変更する（ケース名を付ける）
  1. ケースディレクトリ内の条件ファイル(condition.nml, initial_position.csv)を編集 
  1. `./build/bin/main`で実行。ケース名を入力して計算開始。


## 条件ファイル(condition.nml, initial_position.csv)解説
### condition.nml
  - **リスタート位置 num_restart**
    - 通常は`0`を指定
    - `1以上`にすると、その値に対応するbackupファイルが読み込まれ、そこからリスタートが始まる
  - **初期分布ファイル名 initialDistributionFName**
    - 指定したbackupファイル(.bu)が読み込まれ、それを飛沫初期分布とする
    - 初期分布を固定したくない場合はコメントアウトすること
  - **飛沫周期発生 periodicGeneration**
    - 1秒当たりの発生飛沫数（整数）を指定
    - 初期配置飛沫をすべてNonActiveにしたのち、順次Activateしていくので、初期配置数が飛沫数の上限となる
  - **気流データファイル名 path2FlowFile**
    - 実行ディレクトリからの相対パス、もしくは絶対パスを指定
    - 現在可能な流れ場ファイル：
      - VTK
      - INP
      - FLD
    - CUBE格子(PLOT3D)は、予め非構造格子に変換してから計算してください。
    - .arrayファイルを指定する場合、別途メッシュファイルが必要なので、`meshFile = ***`と指定する
  - **ステップ数オフセット OFFSET**
    - 飛沫計算を、流体連番ファイルの途中の番号から始めたいときに指定
  - **気流データを周期的に用いる場合の先頭と末尾 LoopHead, LoopTail**
    - 任意の区間の流体連番ファイルを繰り返し用いるときに指定（例えば呼吸のサイクル）
    - `(先頭) = (末尾)` とすると、そのステップ数到達後は流れ場の更新が起こらなくなる
    - `(先頭) > (末尾)` とすれば、特殊な処理は起こらず、流体連番ファイルが順番に読み込まれる
### initial_position.csv
  - 初期飛沫の配置帯（直方体）を設定する
  - 左から順に、直方体の中心座標(x,y,z), 直方体の幅(x,y,z)
  - 改行すれば配置帯を複数設定できる


## 方程式

### 飛沫の蒸発方程式

  $$ \frac{dr}{dt} \space = \space -\left(1-\frac{RH}{100}\right) \cdot \frac{D e_{s}(T)}{\rho_{w} R_{v} T} \cdot \frac{1}{r} $$
  
  プログラム内では、２次精度ルンゲクッタ法で解いている。
  
### 飛沫の運動方程式

$$ m \frac{d \mathbf{v}}{dt} \space = \space m \mathbf{g} \space + \space C_D (\mathbf{v}) \space \cdot \space \frac{1}{2} \rho_a S \left | \mathbf{u}_a - \mathbf{v} \right | (\mathbf{u}_a - \mathbf{v}) $$

  プログラム内では、上式を無次元化・離散化した次式を解いている。
    
$$ \bar{\mathbf{v}}^{n + 1} \space = \space \frac{\bar{\mathbf{v}}^{n} \space + \space (\bar{\mathbf{g}} \space + \space C \bar{\mathbf{u}}_a)\Delta \bar{t}}{1 \space + \space C\Delta \bar{t}} \quad \left ( C \space = \space \frac{3 \rho_a}{8 \rho_w} \frac{C_D ( \mathbf{v}^{n} ) \left | \bar{\mathbf{u}}_a - \bar{\mathbf{v}}^{n} \right |}{\bar{r}^{n+1}} \right ) $$

## サブプログラム
  - CUBE2USG
    - CUBE格子を、非構造格子に変換できる
  - droplet2CSV
    - 飛沫計算結果を再度読み込み、統計データ（浮遊数推移など）をCSVファイルに書き出す
  - dropletCount
    - 飛沫計算結果を再度読み込み、カウントボックスを通過した飛沫数を調べる。optionディレクトリ内の"boxList.csv"を、ケースディレクトリに配置する必要がある。
  - initialTranslate
    - 飛沫の初期配置データを読み込み、任意の座標への回転、平行移動を行う。**by Konishi**

## おまけ機能
<<<<<<< HEAD

### **複数ケース連続実行**
実行時にTXTファイル名を入力すると、そのファイルに列挙された複数ケースを連続実行できる
    
### `option/basicSetting.nml`
- adhesionSwitch
  - 飛沫付着判定のスイッチ。OFFにすることないね。要らないね。
- coalescenceLimit
  - 合体が起こらないステップが、指定ステップ数（デフォルトは`10000`）連続した場合、合体判定そのものをOFFにする
  - 合体が起こらないのにいつまでも合体判定し続けるのは計算コストなので、軽減すべく追加した機能
  - `0`を指定すると合体判定オフ
- num_divide
  - 合体判定の分割数（デフォルトは`4`）
  - 全飛沫と合体判定するのは計算コストなので、エリアを分割してそれぞれ別々に合体判定を行う機能
- radiusDistributionFNAME
  - 飛沫初期半径データのファイル名を指定する（デフォルトは咳用のTXTファイル）
  - 咳と会話とでは初期半径が異なるらしく、その場合は別ファイルを指定
  - `data`ディレクトリ内にあることを想定しており、パスではなく単にファイル名を指定
      
### RandomMaker
ランダムな数列を生成するPythonスクリプト。
シミュレーション実行時に毎回乱数を生成するのは比較実験に向かないので、別スクリプトで必要に応じて行う。
- `deadlineMaker.py`
  - ウイルス飛沫の寿命を乱数から生成するスクリプト
- `radiusDistributionMaker.py`
  - 飛沫の初期半径を乱数から生成するスクリプト
      
=======
  - **複数ケース連続実行**
    - 実行時にTXTファイル名を入力すると、そのファイルに列挙された複数ケースを連続実行できる
  - **basicSetting.nml**
    - optionディレクトリ内にある。付着判定のオンオフや、飛沫間合体の設定が可能。初期半径分布ファイルの指定も可能。
    
## CTest
  - コンパイル後、`$ ctest`でCTestが実行可能（buildディレクトリにて）
  - CTestの実行ディレクトリは、`test/`になる（buildディレクトリではない）
  - テスト用プログラムはすべて`test/`で管理しよう
>>>>>>> 800e1971ab14a40e11662faec025ea1b9b7dc95b
