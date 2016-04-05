# If the text is garbled, please convert encoding to UTF-8.
四面体メッシュのFEM弾性体シミュレーション、Spatial Hashingの衝突検出、反応

実行ファイル：Release/FEM.exe

*実行環境等を含めた実行方法:

動作確認：Windows 10　Visual Studio 2013

開発時間と人数: 一年前 一か月　一人

開発言語：C++  OpenGL

外部ライブラリ：glfw glew OpenMesh eigen

実行方法："OpenGL.sln"　Visual Studio 2013で開けたら、F5を押す。

操作方法：

左マウスクリック+ドラッグ　　　視点を回転

スペースボタン　　　　　　　位置をリセット　　　　　　　　　　


*プログラムを作成する上で苦労した箇所は？

デバッグ、FEM弾性体シミュレーションの勉強

*力をいれて作った部分で、「プログラム上」で特に注意してみてもらいたい箇所は？

deformableobject.cppで Spatial HashingのfirstPass（ハッシュの構築）関数、secondPass（検索）関数

*参考にしたソースファイルがあるなら、どの様なところを参考にしましたか？またその部分のファイル名を書いてください

OpenCloth：https://github.com/mmmovania/opencloth
FEM弾性体シミュレーションに参考。

*他人のコードと関数：

deform_other.cpp　弾性体シミュレーションの関連


# FEM_elasticity
This is a Corotation FEM elasticity, Spatial Hashing for collision detetion 
and collision respon implement.

Platform:           Windows 10

Develop Tool:       Visual Studio 2013

Program language:   C++  OpenGL

External libraries: glfw glew OpenMesh eigen

usage:

1.Download zip or clone this repository.

2.Open "OpenGL.sln" using Visual Studio 2013(Only support version 2013 now).

3.Change the Configuration to Release.

4.Run.
In this step you may see "something.dll are not in your compute" something like this,
go "lib" sub folder and copy the "something.dll" file to "Release" folder.After this,
it should work.
  
For Corotation FEM I use the code form https://code.google.com/archive/p/opencloth/.

Screenshots:

![image](https://github.com/duoshengyu/FEM_elasticity/blob/master/screenshots/1.png)

![image](https://github.com/duoshengyu/FEM_elasticity/blob/master/screenshots/2.png)
