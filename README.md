# MESH_CNT

1. This version runs CNT diffusive wave equation.
2. The DT can be variavble.
3. DT is tested at Goodwin as 60,90,120,240,360,600s and they all gives consistant results.
4. This version checks the ratio of Sf between two nodes and applies normal depth incase the ratio is more than 1e5.
5. This version also takes lateral flow at its certain location and applies them properly.