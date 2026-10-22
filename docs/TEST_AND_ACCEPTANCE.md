# 检查与证据

| ID | 检查 | 结果 | 位置 | 限制 |
| --- | --- | --- | --- | --- |
| EVD-0001 | 外部原始文献与官方资料 | STATIC_CHECKED | research/group_meeting_20261022/research/public_sources.json | 43项；部分仅摘要/元数据 |
| EVD-0002 | 标量波动图数值检查 | PASS | research/group_meeting_20261022/evidence/wave_figure_metrics.json | 固定模型SYNTHETIC，无实机 |
| EVD-0003 | 双点反例与多look计算 | PASS | research/group_meeting_20261022/evidence/coherence_metrics.json | 固定seed；独立微态由程序生成 |
| EVD-0004 | 内部交叉审查与修订 | PASS | 受限包scientific_review_*.md | 不是外部同行评审 |
| EVD-0005 | 演示渲染/结构/逐页检查 | PASS | 受限包presentation_validation.json | 版面与文件质量，不是科学性能认证 |
| EVD-0006 | 阅读器数学/链接/离线检查 | PASS | 受限包reader_validation.json | 本地浏览器验证 |
| EVD-0007 | GitHub内容提交与读回 | PASS | docs/turns/receipts/TURN-0001.json；内容提交28f91949b6a2a3b69602cf1630183537eb56245a | 33项提交集合一致，5个关键文件逐字读回一致；外部文件不在仓库 |
| EVD-0008 | 本装置原始栈/复焦扫/独立点靶 | BLOCKED | 尚未提供 | 不得声称核精度/DOF/最高水平通过 |

本轮没有运行主仓库测试，因为不改仿真实现；历史日志中的CI结果仅按历史来源记录，不冒充本轮执行。研究代码有数值收敛、归一化及解析关系对照。实机的预设建议门槛留在完整实验方案，待预实验冻结，尚非验收标准。
