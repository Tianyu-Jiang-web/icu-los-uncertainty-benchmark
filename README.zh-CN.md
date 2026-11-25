# icu-los-uncertainty-benchmark

`drf-replication-ehr-v3`（`scripts.R`、`script_recovery.R`、`39quantiles.R`、`lancet_style.R`、`crps_guidance_lancet.R`）的模块化版本。
功能、随机种子、数据划分、超参网格、分位数网格均与原脚本一致；原文件夹未被修改。

## 目录结构

```
config.R                    路径 / 分位数网格 / 开关（可用环境变量覆盖）
run_train.R                 流程1：完整训练（原 scripts.R）
run_recovery.R              流程2：torch 模型恢复（原 script_recovery.R）
run_predict_visualize.R     流程3：39 分位预测 + 全部图表（原 39quantiles.R + crps_guidance_lancet.R）
R/
  pipeline.R                load_modules()、select_torch_device()、run_stages()
  utils_metrics.R           纯函数：pinball、CRPS、加权分位数 wquant_matrix
  packages.R                三条流程各自的 library() 顺序
  data_prep.R               特征构建 + 时间切分（三条流程共用）
  eda_plots.R / table1.R    EDA 图 / Table 1
  training_setup.R          CV 折与 19 分位网格
  models/                   每个模型一个文件：drf qrf ranger xgboost ens mcd ddnn ensemble_drf bnn hybrid + torch_common
                            （每个文件含 stage_train_* 与 stage_predict39_*）
  predict_setup.R           39 分位网格、加载 rds/.pt、汇总 all_q / all_qgrid
  analysis/                 evaluation heterogeneity fairness calibration sensitivity
  plots/                    heterogeneity_plots fairness_plots sensitivity_plots crps_guide
  outputs.R                 保存超参 / rds / 汇总 csv / 包版本
  recovery.R                流程2 的各阶段
  plot_lancet_style.R       Lancet 风格绘图工具（原 lancet_style.R）
```

## 运行

在本目录下：

```bash
Rscript run_predict_visualize.R     # 用已保存模型重新预测并生成全部图表（默认读取原项目 figures/final）
Rscript run_recovery.R              # 重训 5 个 torch 模型并保存 .pt
Rscript run_train.R                 # 完整调参 + 训练（耗时数小时）
```

每个 `stage_*` 函数就是原脚本的一段代码（逐字保留），按顺序在全局环境里执行，因此数据在阶段之间的流动与原脚本相同。
可只运行部分阶段：

```bash
DRF_ONLY="stage_prepare_features,stage_split_data" Rscript run_predict_visualize.R
DRF_SKIP="stage_eval39_sensitivity,stage_plot39_sensitivity" Rscript run_predict_visualize.R
```

## 路径与环境变量

| 变量 | 默认值 | 说明 |
|---|---|---|
| `DRF_DATA_PATH` | `../drf-replication-ehr-v3/data/FINAL_master_dataset_9Feb2026(in).csv` | 主数据 |
| `DRF_OUTPUT_DIR` | `./outputs` | 本项目输出目录 |
| `DRF_MODELS_DIR` | `outputs/final`（若已有 rds），否则 `../drf-replication-ehr-v3/figures/final` | `all_model_results.rds` 所在位置 |
| `DRF_LEGACY_DIR` | `../drf-replication-ehr-v3` | 原项目 |
| `DRF_FORCE_CPU=1` | 未设置 | 不使用 Apple `mps`，用于可比对的测试 |

## 与原脚本的差异（仅此几处）

1. 硬编码的绝对路径 / `data/...` 相对路径改为读取 `cfg`。
2. 设备选择 `torch_device(if (mps) "mps" else "cpu")` 改为 `select_torch_device()`（默认行为相同）。
3. 三个脚本里重复的数据预处理、`crps_from_quantiles`、`wquant_matrix` 合并为一份（已验证与三个原版结果逐位相同）。
4. `crps_guidance_lancet.R` 原本需手动在 `39quantiles.R` 之后运行，现作为 `run_predict_visualize.R` 的最后一个阶段（`cfg$run_crps_guide`）。
5. 阶段函数内的代码不会在控制台自动打印未赋值表达式（不影响任何结果与输出文件）。

训练流程（`scripts.R`）与恢复流程（`script_recovery.R`）里各自定义的神经网络函数默认参数略有不同（如默认 epochs），为保证结果一致，保留各自版本，没有合并。
