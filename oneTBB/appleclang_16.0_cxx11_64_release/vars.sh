#!/bin/sh
#
# Copyright (c) 2005-2021 Intel Corporation
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

export TBBROOT=/Users/zhangwenqian/UNSW/pivoter/oneTBB

DYLD_LIBRARY_PATH="/Users/zhangwenqian/UNSW/pivoter/oneTBB/appleclang_16.0_cxx11_64_release:${DYLD_LIBRARY_PATH}"; export DYLD_LIBRARY_PATH
LIBRARY_PATH="/Users/zhangwenqian/UNSW/pivoter/oneTBB/appleclang_16.0_cxx11_64_release:${LIBRARY_PATH}"; export LIBRARY_PATH
CPATH="${TBBROOT}/include:${CPATH}"; export CPATH
PKG_CONFIG_PATH="/Users/zhangwenqian/UNSW/pivoter/oneTBB/appleclang_16.0_cxx11_64_release/pkgconfig:${PKG_CONFIG_PATH}"; export PKG_CONFIG_PATH


