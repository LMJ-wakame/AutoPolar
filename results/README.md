# 结果文件说明

根目录的 algorithm_derivation.json、spc_study.json、gpc2_study.json、gpc4_study.json、gpc4_cost_audit.json 是迁移前已经保存的历史结果，保留原值。

cpp_validation/ 保存2026-10-25本次C++运行的结果与测试日志。新旧运行的随机引擎不同，不应按相同随机样本逐项解释其平均操作量差异。验证的目标是各实现相对于精确基线的一致性。

未加权操作数不等于运行时间；FER演示中的少量帧不是正式大规模性能基准。
