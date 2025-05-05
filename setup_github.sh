#!/bin/bash
# 初始化git存儲庫並準備上傳到GitHub

# 設置您的GitHub用戶名和存儲庫名稱（請替換為您的信息）
GITHUB_USERNAME="your-username"
REPO_NAME="multi-tumor-network"

# 初始化git存儲庫
git init

# 添加文件
git add multi_tumor_network.py
git add README.md
git add LICENSE
git add .gitignore
git add requirements.txt
git add USAGE.md
git add setup_github.bat
git add setup_github.sh
git add docs/*.md

# 創建output目錄（如果不存在）
mkdir -p output

# 添加數據文件（如果有）
git add -f data/*.csv

# 排除output目錄
git rm -r --cached output 2>/dev/null || true

# 提交更改
git commit -m "Initial commit: Multi-Tumor Network Visualization Tool"

# 添加遠程存儲庫
git remote add origin https://github.com/$GITHUB_USERNAME/$REPO_NAME.git

# 推送到GitHub
echo "準備好推送到GitHub。運行以下命令完成上傳："
echo "git push -u origin master"

echo ""
echo "如果您尚未創建GitHub存儲庫，請先訪問 https://github.com/new 創建一個名為 $REPO_NAME 的存儲庫" 