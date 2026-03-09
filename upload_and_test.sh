#!/bin/bash

# 一键上传和测试脚本 - 在本地运行

echo "=========================================="
echo "Uploading and Testing on Server"
echo "=========================================="
echo ""

SERVER="tods2"
REMOTE_DIR="~/pivoter"

echo "Target server: $SERVER"
echo "Remote directory: $REMOTE_DIR"
echo ""

# 1. 上传代码
echo "[1/3] Uploading code to server..."
scp -r /Users/zhangwenqian/UNSW/pivoter $SERVER:~/ 2>&1 | grep -v "Warning"
if [ $? -eq 0 ]; then
    echo "✓ Upload successful"
else
    echo "✗ Upload failed"
    exit 1
fi
echo ""

# 2. 在服务器上运行部署脚本
echo "[2/3] Running deployment on server..."
echo "This may take 5-10 minutes..."
echo ""

ssh $SERVER << 'ENDSSH'
cd ~/pivoter
chmod +x deploy_server.sh test_server_32threads.sh
./deploy_server.sh
ENDSSH

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Server testing complete"
else
    echo ""
    echo "✗ Server testing failed"
    exit 1
fi
echo ""

# 3. 下载结果
echo "[3/3] Downloading results..."
scp $SERVER:~/pivoter/server_performance_results_*.txt . 2>&1 | grep -v "Warning"
if [ $? -eq 0 ]; then
    echo "✓ Results downloaded"
    echo ""
    echo "Results file:"
    ls -lh server_performance_results_*.txt | tail -1
else
    echo "Note: Could not download results (may not exist yet)"
fi

echo ""
echo "=========================================="
echo "Complete!"
echo "=========================================="
echo ""
echo "To view results on server:"
echo "  ssh $SERVER"
echo "  cat ~/pivoter/server_performance_results_*.txt"
