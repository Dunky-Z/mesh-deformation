
#### ASAP
在三维变形操作中， 为了满足特定限制，如控制点位置，那么拉伸和扭曲就不可避免。因此在三维模型变形操作如果能在局部发生尽可能少的拉伸和扭曲，那么细节就可以被保持。这种变形算法称为ARAP（As Rigid As Possible）算法。
###### 先只考虑局部情况
让$c_i$表示对应于顶点$i$的局部网格，也就是此顶点相邻的所有面。$c_i$变形后表示为$c'_i$。假如变形是rigid的，也就是没有拉伸扭曲，那么有：

$$
p'_i-p'_j=\bm{R_i}(p_i-p_j),\forall j\in N(i)
$$

若变形不是刚性的那么可以定义以下能量函数，求得$R$使得能量函数值最小：

$$
E(c_i,c'_i)=\sum_{j \in N(i)} w_{ij} \Vert (p'_i-p'_j)-\bm{R_i}(p_i-p_j) \Vert^2
$$

为简洁设$e_{ij} = p_i-p_j$,将能量函数用矩阵形式改写：

$$
\sum_{j \in N(i)} w_{ij} (e^{'T}_{ij}e'_{ij}-2e^{'T}_{ij}\bm{R_i}e_{ij}+e^{T}_{ij})
$$

为求最小值，只需要最小化中间含$\bm{R_i}$的项，即

$$
arg\min_{R_i}\sum_j -2w_{ij}e^{'T}_{ij}R_ie_{ij}=argmax_{R_i}Tr(\bm{R_i}\sum_j w_{ij}e_{ij}e^{'T}_{ij})
$$

由于当*$M$是个**对称半正定矩阵**时，对于任意正交矩阵$R$,都有$Tr(M)\geq Tr(RM)$*。故上式中$\bm{R_i}\sum_j w_{ij}e_{ij}e^{'T}_{ij}$为半正定时，其$Tr$最大。此时可由矩阵$\sum_j w_{ij}e_{ij}e^{'T}_{ij}$奇异值分解求出$\bm{R_i}$：

$$
\sum_j w_{ij}e_{ij}e^{'T}_{ij},\quad \bm{R_i}=V_iU^T_i
$$

##### 整体的变形
这个网格上的能量函数可以由单个网格上能量函数之和得到：

$$
E\left(\mathcal{S}^{\prime}\right)=\sum_{i=1}^{n}w_{i}E\left(\mathcal{C}_{i},\mathcal{C}_{i}^{\prime}\right)=\\=\sum_{i=1}^{n}w_{i}\sum_{j\in\mathcal{N}(i)}w_{ij}\left\|\left(\mathbf{p}_{i}^{\prime}-\mathbf{p}_{j}^{\prime}\right)-\mathbf{R}_{i}\left(\mathbf{p}_{i}-\mathbf{p}_{j}\right)\right\|^{2}
$$

能量函数的最小值可以通过对变形后顶点的位置求导得到：

$$
\frac{\partial{E}\left(\mathcal{S}^{\prime}\right)}{\partial\mathbf{p}_{i}^{\prime}}=\frac{\partial}{\partial\mathbf{p}_{i}^{\prime}}(\sum_{j\in\mathcal{N}(i)}w_{ij}\left\|\left(\mathbf{p}_{i}^{\prime}-\mathbf{p}_{j}^{\prime}\right)-\mathbf{R}_{i}\left(\mathbf{p}_{i}-\mathbf{p}_{j}\right)\right\|^{2}+\\\sum_{j\in\mathcal{N}(i)}w_{ji}\left\|\left(\mathbf{p}_{j}^{\prime}-\mathbf{p}_{i}^{\prime}\right)-\mathbf{R}_{j}\left(\mathbf{p}_{j}-\mathbf{p}_{i}\right)\right\|^{2})
\\=\sum_{j\in\mathcal{N}(i)}2w_{ij}\left(\left(\mathbf{p}_{i}^{\prime}-\mathbf{p}_{j}^{\prime}\right)-\mathbf{R}_{i}\left(\mathbf{p}_{i}-\mathbf{p}_{j}\right)\right)+\\\sum_{j\in\mathcal{N}(i)}-2w_{ji}\left(\left(\mathbf{p}_{j}^{\prime}-\mathbf{p}_{i}^{\prime}\right)-\mathbf{R}_{j}\left(\mathbf{p}_{j}-\mathbf{p}_{i}\right)\right)
$$
由于$w_{ij}=w_{ji}$，可得
$$
\frac{\partial{E}\left(\mathcal{S}^{\prime}\right)}{\partial\mathbf{p}_{i}^{\prime}}=\sum_{j\in\mathcal{N}(i)}4w_{ij}\left(\left(\mathbf{p}_{i}^{\prime}-\mathbf{p}_{j}^{\prime}\right)-\frac{1}{2}\left(\mathbf{R}_{i}+\mathbf{R}_{j}\right)\left(\mathbf{p}_{i}-\mathbf{p}_{j}\right)\right)
$$
上式等于0，可得：
$$
\sum_{j\in\mathcal{N}(i)}w_{ij}\left(\mathbf{p}_{i}^{\prime}-\mathbf{p}_{j}^{\prime}\right)=\sum_{j\in\mathcal{N}(i)}\frac{w_{ij}}{2}\left(\mathbf{R}_{i}+\mathbf{R}_{j}\right)\left(\mathbf{p}_{i}-\mathbf{p}_{j}\right)
$$
化成矩阵形式：
$$
\bm{Lp'}=\bm{b}
$$

#### 非线性系统局部最优问题
上述系统中，未知数是$p',R$。但是$R$又需要$p'$才能求得。自然我们想到一种迭代方法：我们固定$p',R$其中一组变量。优化出使得能量函数最小的另一组变量的值，然后将求得的这组值再固定，进入下一次迭代去优化之前固定的那组变量，如此往复，直至最终收敛。在这里先算$p'$还是$R$也是个问题。
