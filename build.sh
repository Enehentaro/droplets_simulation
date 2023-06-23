#!/bin/bash

#ビルド用シェルスクリプト

mkdir -v build
cd build
cmake .. -D CMAKE_BUILD_TYPE=debug
make
#ctest
cd ..


#実行ファイルへのパスが、環境変数パスになければ追加
#あくまで子プロセスの環境変数を書き換えるにすぎない
# $ source build.sh で実行すると、子プロセス終了後も環境変数は引き継がれる
path2bin=$(pwd)/build/bin
if [[ ! "$PATH" =~ $path2bin ]]; then
    export PATH=$PATH:$path2bin
fi
