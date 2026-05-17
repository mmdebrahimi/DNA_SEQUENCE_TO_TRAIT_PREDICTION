@echo off
cd /d C:\Users\Farshad\PythonProjects\dna_decode
set HF_HOME=D:\hf_cache
uv run python -u scripts\stage1_n40_cipro.py > stage1b_detached.log 2>&1
