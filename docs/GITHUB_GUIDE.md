# GitHub上傳指南

本文檔提供如何將多腫瘤網絡視覺化工具上傳到GitHub的詳細步驟。

## 前置條件

1. 安裝Git工具：https://git-scm.com/downloads
2. 創建GitHub帳號：https://github.com/join
3. 配置Git身份：
   ```bash
   git config --global user.name "您的名字"
   git config --global user.email "您的電子郵件"
   ```

## 在GitHub創建新存儲庫

1. 登錄GitHub，訪問：https://github.com/new
2. 輸入存儲庫名稱：`multi-tumor-network`
3. 添加描述（可選）：`多腫瘤網絡視覺化工具`
4. 選擇公共(Public)或私有(Private)存儲庫
5. 不要初始化存儲庫（不選擇添加README、.gitignore或許可證）
6. 點擊「創建存儲庫」

## 使用提供的腳本上傳

### Windows用戶

1. 編輯`setup_github.bat`文件，將`your-username`替換為您的GitHub用戶名
2. 運行腳本:
   ```
   setup_github.bat
   ```
3. 根據提示運行:
   ```
   git push -u origin master
   ```
4. 輸入您的GitHub用戶名和密碼或令牌（如果已配置SSH密鑰則不需要）

### macOS/Linux用戶

1. 編輯`setup_github.sh`文件，將`your-username`替換為您的GitHub用戶名
2. 添加執行權限：
   ```bash
   chmod +x setup_github.sh
   ```
3. 運行腳本:
   ```bash
   ./setup_github.sh
   ```
4. 根據提示運行:
   ```bash
   git push -u origin master
   ```
5. 輸入您的GitHub用戶名和密碼或令牌（如果已配置SSH密鑰則不需要）

## 手動上傳步驟

如果您不想使用腳本，可以按照以下步驟手動操作：

1. 初始化git存儲庫：
   ```bash
   git init
   ```

2. 添加文件：
   ```bash
   git add multi_tumor_network.py README.md LICENSE .gitignore requirements.txt USAGE.md
   git add -f data/*.csv  # 添加數據文件（如果有）
   ```

3. 提交更改：
   ```bash
   git commit -m "Initial commit: Multi-Tumor Network Visualization Tool"
   ```

4. 添加遠程存儲庫：
   ```bash
   git remote add origin https://github.com/您的用戶名/multi-tumor-network.git
   ```

5. 推送到GitHub：
   ```bash
   git push -u origin master
   ```

## 更新現有存儲庫

如果您需要更新已上傳的存儲庫：

1. 添加修改的文件：
   ```bash
   git add 修改的文件
   ```

2. 提交更改：
   ```bash
   git commit -m "更新說明信息"
   ```

3. 推送到GitHub：
   ```bash
   git push
   ```

## 常見問題

**Q: 上傳時出現「Permission denied」錯誤怎麼辦？**  
A: 確認您的GitHub帳號有權訪問該存儲庫，或者檢查您的身份驗證信息。

**Q: 如何使用個人訪問令牌代替密碼？**  
A: 訪問GitHub設置 -> Developer settings -> Personal access tokens生成令牌，然後在推送時使用該令牌作為密碼。

**Q: 如何忽略某些文件不上傳？**  
A: 編輯`.gitignore`文件，添加您想忽略的文件或目錄。

**Q: 推送時發生衝突怎麼辦？**  
A: 先執行`git pull`拉取最新更改，解決衝突後再推送。

**Q: 如何查看存儲庫狀態？**  
A: 使用命令`git status`查看當前存儲庫狀態。 