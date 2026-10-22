# 任务DAG

写入边界：研究输出与治理；不修改cop_oct_sim或既有测试。DONE只表示该项声明的研究交付完成。

| ID | 工作 | 依赖 | 负责人 | 输出 | 证据 | 状态 |
| --- | --- | --- | --- | --- | --- | --- |
| TASK-0001 | 来源核验与文献分类 | 无 | 研究代理 | public_sources、references | EVD-0001 | DONE |
| TASK-0002 | 主讲与条件推导交叉检查 | TASK-0001 | 根与研究代理 | 完整文本；公开治理摘要 | EVD-0004 | DONE |
| TASK-0003 | 科学仿真与图表 | TASK-0001 | 波动/OCT代理 | 公开code与metrics | EVD-0002/0003 | DONE |
| TASK-0004 | 42页可编辑演示与PDF | TASK-0002/0003 | 演示代理 | 外部演示文件 | EVD-0005 | DONE |
| TASK-0005 | 离线阅读器与文件夹 | TASK-0002/0003 | 阅读器代理与根 | 外部完整包 | EVD-0006 | DONE |
| TASK-0006 | 原子同步与读取核实 | TASK-0004/0005 | 根 | 治理+公共子集 | EVD-0007 | DONE_PUBLIC_SUBSET |
| TASK-0007 | 独立实机核与DOF验收 | TASK-0002 | 后续实验 | 尚无原始数据 | EVD-0008 | EXTERNAL_EVIDENCE_MISSING |
