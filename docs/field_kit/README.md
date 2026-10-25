# 现场工具包 / Field kit

路线 2 需要真实楼梯。物理采集只能由你（或你的同学/合作者）完成；这个工具包把现场工作压到最低。

## 包含

| 文件 | 用途 |
|---|---|
| `permission_letter_zh.md` | 给建筑管理/文保部门的许可申请信模板 |
| `count_sheet_A4.pdf` | 14 天人流计数表，A4 打印 |
| `scale_bar_A4.pdf` | 100 mm 比例尺 + 50 mm 棋盘格 + 编码点，A4 打印 |
| `photogrammetry_guide.md` | 手机拍摄、导出、命名、常见错误 |
| `../../../data/field/count_template.csv` | 计数数据录入模板 |

## 最小工作量

- **3 座标定楼梯**：优先选有门禁/闸机计数的图书馆、办公楼、教学楼；能导出计数就不用手工数 14 天。
- **每座 5 级踏步**：手机拍 40–80 张/级，10–15 分钟/楼梯。
- **许可**：用模板改一页，盖章/签名。
- **总计**：每座楼梯现场 1–2 小时（含许可、拍摄、记录）；如果门禁计数可导出，14 天计数可以省掉。

## 交付给我

```
q1_rebuild/data/field/raw/
├── C1.zip          # 3 座标定楼梯
├── C2.zip
├── C3.zip
├── A1.zip          # 3 座应用楼梯（历史建筑）
├── A2.zip
├── A3.zip
├── metadata.csv    # 按 data/field/README.md 的 schema
└── counts.csv      # 按 count_template.csv
```

## 我做完剩下的

- 三维重建检查/剖面提取（`mesh_slice.py`）
- 参考面拟合、磨深计算
- 标定 κ、反演 N/δ/车道数、置信区间
- 真实案例图、Results 小节、摘要/cover letter 更新
- 重编译 PDF/DOCX

## 如果拿不到 3 座

2 座也能做，区间会变宽；1 座只能做 demo，不能做标定。不要用编的数据凑数。
