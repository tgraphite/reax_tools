#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ReaxTools Reporter - 将test_data文件内容嵌入到HTML模板中
生成包含内嵌数据的report.html文件，并内嵌所有CSS和JS文件
"""

import os
import json
import re
import requests
from pathlib import Path
from urllib.parse import urlparse


def read_file_content(file_path):
    """读取文件内容"""
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            return f.read()
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None


def download_url_content(url):
    """下载URL内容"""
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        return response.text
    except Exception as e:
        print(f"Error downloading {url}: {e}")
        return None


def escape_js_string(content):
    """转义JavaScript字符串中的特殊字符"""
    # 转义JavaScript字符串中的特殊字符
    content = content.replace("\\", "\\\\")
    content = content.replace('"', '\\"')
    content = content.replace("\n", "\\n")
    content = content.replace("\r", "\\r")
    content = content.replace("\t", "\\t")
    return content


def create_data_json(test_data_dir):
    """创建包含所有数据的JSON对象"""
    data_files = [
        "species_count.csv",
        "bond_count.csv",
        "atom_bonded_num_count.csv",
        "ring_count.csv",
        "reactions.dot",
        "key_molecules_reactions.csv",
    ]

    data_dict = {}

    for filename in data_files:
        file_path = test_data_dir / filename
        if file_path.exists():
            content = read_file_content(file_path)
            if content is not None:
                data_dict[filename] = content
                print(f"✓ Loaded {filename}")
            else:
                print(f"✗ Failed to load {filename}")
        else:
            print(f"✗ File not found: {filename}")

    return data_dict


def extract_and_inline_resources(template_content, template_dir):
    """提取并内嵌所有CSS和JS资源"""
    print("Processing external resources...")

    # 处理CSS文件
    css_pattern = r'<link\s+rel="stylesheet"\s+href="([^"]+)"'
    css_matches = re.findall(css_pattern, template_content)

    for css_href in css_matches:
        print(f"Processing CSS: {css_href}")
        css_content = None

        if css_href.startswith("http"):
            # CDN文件
            css_content = download_url_content(css_href)
        else:
            # 本地文件
            css_path = template_dir / css_href.lstrip("./")
            if css_path.exists():
                css_content = read_file_content(css_path)

        if css_content:
            # 替换link标签为内嵌style标签
            old_tag = f'<link rel="stylesheet" href="{css_href}">'
            new_tag = f"<style>\n{css_content}\n</style>"
            template_content = template_content.replace(old_tag, new_tag)
            print(f"✓ Inlined CSS: {css_href}")
        else:
            print(f"✗ Failed to inline CSS: {css_href}")

    # 处理JavaScript文件
    js_pattern = r'<script\s+src="([^"]+)"[^>]*></script>'
    js_matches = re.findall(js_pattern, template_content)

    for js_src in js_matches:
        print(f"Processing JavaScript: {js_src}")
        js_content = None

        if js_src.startswith("http"):
            # CDN文件
            js_content = download_url_content(js_src)
        else:
            # 本地文件
            js_path = template_dir / js_src.lstrip("./")
            if js_path.exists():
                js_content = read_file_content(js_path)

        if js_content:
            # 替换script标签为内嵌script标签
            old_tag = f'<script src="{js_src}"></script>'
            new_tag = f"<script>\n{js_content}\n</script>"
            template_content = template_content.replace(old_tag, new_tag)
            print(f"✓ Inlined JavaScript: {js_src}")
        else:
            print(f"✗ Failed to inline JavaScript: {js_src}")

    return template_content


def generate_report_html(template_path, output_path, data_dict, project_id="未定义"):
    """生成包含内嵌数据的HTML报告"""
    try:
        # 读取模板文件
        with open(template_path, "r", encoding="utf-8") as f:
            template_content = f.read()

        # 内嵌所有CSS和JS资源
        template_dir = template_path.parent
        template_content = extract_and_inline_resources(template_content, template_dir)

        # 将数据转换为JSON字符串
        data_json = json.dumps(data_dict, ensure_ascii=False, indent=2)

        # 替换占位符
        placeholder = "<!-- DATA_PLACEHOLDER -->"
        if placeholder in template_content:
            template_content = template_content.replace(placeholder, data_json)
            print(f"✓ Data placeholder replaced")
        else:
            print(f"✗ Data placeholder not found in template")
            return False

        project_id_placeholder = "DATA_PROJECT_ID"
        if project_id_placeholder in template_content:
            template_content = template_content.replace(
                project_id_placeholder, project_id
            )
            print(f"✓ Project ID placeholder replaced")
        else:
            print(f"✗ Project ID placeholder not found in template")
            return False

        # 写入输出文件
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(template_content)

        print(f"✓ Report generated: {output_path}")
        return True

    except Exception as e:
        print(f"✗ Error generating report: {e}")
        return False


def main():
    """主函数"""
    # 改为生产环境
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--input_dir", type=str, help="测试数据目录")
    parser.add_argument(
        "-p", "--project_id", type=str, help="项目ID", default="您尚未设置项目ID"
    )
    args = parser.parse_args()
    input_dir = Path(args.input_dir)
    output_path = input_dir / "report.html"

    # 获取当前脚本所在目录
    script_dir = Path(__file__).parent
    template_path = script_dir / "template.html"

    print("ReaxTools Reporter - Enhanced Version")
    print("=" * 60)
    print("Features: Data embedding + CSS/JS inlining")
    print("=" * 60)

    # 检查必要文件是否存在
    if not template_path.exists():
        print(f"✗ Template file not found: {template_path}")
        return

    if not input_dir.exists():
        print(f"✗ Test data directory not found: {input_dir}")
        return

    print(f"Template: {template_path}")
    print(f"Input Data: {input_dir}")
    print(f"Output: {output_path}")
    print("-" * 60)

    # 创建数据字典
    print("Loading data files...")
    data_dict = create_data_json(input_dir)

    if not data_dict:
        print("✗ No data files could be loaded")
        return

    print(f"✓ Loaded {len(data_dict)} data files")
    print("-" * 60)

    # 生成HTML报告
    print("Generating HTML report with inlined resources...")
    if generate_report_html(template_path, output_path, data_dict, args.project_id):
        print("=" * 60)
        print("✓ Report generation completed successfully!")
        print(f"📁 Output file: {output_path}")
        print("🌐 You can now open report.html in your browser")
        print("📊 All data is embedded and no external requests needed")
        print("🎨 All CSS and JavaScript files are inlined")
        print("📦 The report is completely self-contained")
    else:
        print("✗ Report generation failed")


if __name__ == "__main__":
    main()
