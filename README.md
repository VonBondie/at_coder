# 環境構築
## online-judge-toolsのインストール
``` shell(online-judge-toolsのインストール)
pip3 install online-judge-tools

# ojのインストールチェック
oj --version
```

## libxml2及びlibxsltが無いエラーが出たときは
``` shell
sudo apt install libxml2-dev libxslt-dev python3-dev
```

## atcoder-cliのインストール
atcoder-cliはnpmによって提供されているので、npmをインストールするためにnvmをインストールする
[fish 環境にnvm + fish-nvmを導入した時のメモ | DevelopersIO](https://dev.classmethod.jp/articles/fish-nvm/)
``` shell
curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.35.0/install.sh | bash
# fish-nvmのインストール
fisher install jorgebucaran/fish-nvm
nvm -v
```

nodeのインストール
[nvmを使ってNode.jsをインストールする #Node.js - Qiita](https://qiita.com/pyon_kiti_jp/items/da5080e9c7454e935aeb)
[Node.js](https://nodejs.org/en)でバージョンを調べ、以下のコマンドを打つ
``` shell
nvm install [version]
```
atcoder-cliのインストール
``` shell
npm install -g atcoder-cli
```
## atcoder-cliの設定
#### ojのパスを設定
``` bash
acc check-oj
# online-judge-tools is availableと表示されればOK
```
##### テストケースのディレクトリ名の設定
``` bash
cd `acc config-dir` #fishでは対応していなさげなのでよしなにやる

# config.dirがあるので、該当業を以下のように書き換える
"default-test-dirname-format": "test"
```

##### テンプレートファイルの配置
``` bash
cd `acc config-dir`
mkdir cpp
cd cpp

cp なにがしテンプレート .

touch template.json
# template.jsonの中身に以下の内容を記述
{
  "task": {
    "program": ["main.cpp"],
    "submit": "main.cpp"
  }
}

# config.jsonを編集、該当業を以下のように書き換える
"default-template": "cpp"
```

# コマンドメモ

## コンテストディレクトリの作成&テストデータのフェッチ
acc new [contestname] --choice all

## ログイン関連
acc login  // accで再ログインが必要になったとき  
oj login https://atcoder.jp/  // ojでログインが必要になったとき  
