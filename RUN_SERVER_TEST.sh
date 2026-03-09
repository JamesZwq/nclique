#!/bin/bash
# 完整的服务器测试命令 - 直接复制粘贴到终端执行

echo "=========================================="
echo "开始服务器测试"
echo "=========================================="
echo ""

# 上传代码
echo "[1/3] 上传代码到服务器..."
cd /Users/zhangwenqian/UNSW/pivoter
scp -r . tods2:~/pivoter/

echo ""
echo "[2/3] 在服务器上编译和测试..."

# SSH到服务器并执行
ssh tods2 << 'ENDSSH'
cd ~/pivoter
echo "当前目录: $(pwd)"
echo "文件列表:"
ls -lh deploy_server.sh test_server_32threads.sh

echo ""
echo "设置权限..."
chmod +x deploy_server.sh test_server_32threads.sh

echo ""
echo "开始部署和测试..."
./deploy_server.sh

echo ""
echo "测试完成！"
ENDSSH

echo ""
echo "[3/3] 下载结果..."
scp tods2:~/pivoter/server_performance_results_*.txt .

echo ""
echo "=========================================="
echo "完成！查看结果："
echo "=========================================="
ls -lh server_performance_results_*.txt
cat server_performance_results_*.txt | grep -E "(Speedup|Testing with)"

echo ""
echo "完整结果已保存到本地"
