## 任务说明
ReaxTools是一个高性能的C++反应动力学模拟后处理软件，输出各类纯文本分析结果。
编写其对应的可视化组件，基于Python，读取ReaxTools的原始输出文件，根据模版制作成单一HTML文件（嵌入js和css），从而为用户提供本地、可交互、易读的报告。

## 可用资源
file_description.md: 一份ReaxTools输出文件格式的详细说明
web_version_reference/*: ReaxTools有一个网页版（基于WASM），网页版自带基础可视化，你可以参考其中的一部分逻辑

## 需编码的文件
template.html: HTML空白模版，可视化组件的用途就是用实际数据填充它，可以直接导入现成的css和js代码，而不需要复制，以避免内容太多。这个模版本身不必能够直接打开，只是模版。
reax_reporter.py: 用于从数据和模版制作report.html的代码。

## 编码风格
- 简明、含义确凿的注释
- python使用小写下划线命名，javascript使用小写驼峰命名，符合语言的一般习惯
- 可调参数提到最前面，并加上简要注释

## 使用方式
使用python reax_reporter.py [目录名] 可以扫描结果目录下的文件，并在目录中生成报告。