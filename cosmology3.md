# 宇宙学基本方程
宇宙学本质上是广义相对论与统计力学相结合的应用。在这里，唯一有效的长程力是引力；并且我们对单个粒子的行为不感兴趣，而是对大量粒子的集体行为感兴趣。因此，我们需要两个方程：
1. 爱因斯坦场方程，描述引力。
2. 波尔兹曼方程，描述物质与辐射的统计力学。
## 3.1 爱因斯坦场方程
爱因斯坦方程中，该方程将描述几何形状的爱因斯坦张量与物质的能量-动量张量联系起来：
$$G_{\mu\nu} + \Lambda g_{\mu\nu} = 8\pi G T_{\mu\nu}\tag{3.1}$$
$$G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2} R g_{\mu\nu}\tag{3.2}$$
其中 \(G_{\mu\nu}\) 是爱因斯坦张量；\(R_{\mu\nu}\) 是里奇张量，依赖于度规及其导数；\(R\) 是标量曲率，是里奇张量的收缩(\(R = g^{\mu\nu} R_{\mu\nu}\))，\(g_{\mu\nu}\) 是度规张量，\(T_{\mu\nu}\) 是能量-动量张量，\(G\) 是引力常数。一般认为式\((3.1)\)左侧是一个度规的函数，右侧代表宇宙中包含的物质，这样爱因斯坦方程就将二者联系在了一起。这个方程形式上是简单的，但内部蕴含着非常复杂的物理，宇宙及其宇宙内物质的演化都与此相关。

里奇张量可以用克里斯托弗符号表示：
$$R_{\mu\nu} = \Gamma^\alpha_{\mu\nu,\alpha} - \Gamma^\alpha_{\mu\alpha,\nu} + \Gamma^\alpha_{\beta\alpha}\Gamma^\beta_{\mu\nu} - \Gamma^\alpha_{\beta\nu}\Gamma^\beta_{\mu\alpha}\tag{3.3}$$
我们先来考虑FLRW度规下的平直时空：此时由于克氏符正比于\(\dot{a}\)，所以右式前两项正比于\(\ddot{a}\)，后两项正比于\((\dot{a})^2\)。

考虑到R表示空间曲率，是几何不变量，不应该随\(a(t_0)\)的值改变而改变，所以我们令\(a(t_0)=1\)，则右式前两项正比于\(\ddot{a}/a\)，后两项正比于\((\dot{a}/a)^2 = H^2\)。

我们先考虑\(\mu=\nu=0\)，此时
$$R_{00} = \Gamma^\alpha_{00,\alpha} - \Gamma^\alpha_{0\alpha,0} + \Gamma^\alpha_{\beta\alpha} \Gamma^\beta_{00} - \Gamma^\alpha_{\beta0} \Gamma^\beta_{0\alpha}$$

由于当克氏符下标都为0时克氏符为0,所以第一项和第三项为0；又由于\(\alpha\)和\(\beta\)都是空间坐标，所以我们用\(i,j\)来代替，则上式变为：
$$R_{00} = -\Gamma^i_{0i,0} - \Gamma^i_{j0} \Gamma^j_{0i}$$
带入克氏符计算得到：
$$R_{00} = -3\frac{\ddot{a}}{a}\tag{3.4}$$
空间部分同理计算得到：
$$R_{ij} = \delta_{ij} [2\dot{a}^2 + a\ddot{a}]\tag{3.5}$$
由此，里奇标量为：
$$R \equiv g^{\mu\nu} R_{\mu\nu} = -R_{00} + \frac{1}{a^2} R_{ii} = 6[\frac{\ddot{a}}{a} + (\frac{\dot{a}}{a})^2]\tag{3.6}$$
对于式\(3.1\)，如果把宇宙常数移到右边，我们就可以定义宇宙常数对应力能量张量的贡献：
$$T_{(\Lambda)^\mu_\nu} = - \frac{\Lambda}{8\pi G}\delta^\mu_\nu = \begin{pmatrix}-\rho_\Lambda&0&0&0\\0&-\rho_\Lambda&0&0\\0&0&-\rho_\Lambda&0\\0&0&0&-\rho_\Lambda\end{pmatrix}, where \rho_\Lambda \equiv \frac{\Lambda}{8\pi G} \tag{3.7}$$
其中\(\rho_\Lambda\)是宇宙常数的有效能量密度，并且我们可以得到\(P_\Lambda = -\rho_\Lambda\)，即\(\omega_\Lambda = -1\)。在之前的第二章中，我们知道\(\omega_\Lambda\)取其他值时，宇宙常数的有效能量密度都不会保持恒定。所以我们将宇宙常数移到右侧最大的好处就是可以轻松的将爱因斯坦场方程的结果推广到所有非宇宙常数（动态变化的）暗能量上。

继续推倒均匀宇宙中的爱因斯坦场方程，取时间部分：
$$R_{00} - \frac{1}{2}g_{00}R = 8\pi G T_{00}$$
其中\(T_{00}\)可以认为是宇宙中所有物质家在一起的能量密度，即\(\rho\)。在带入相关数值，得到：
$$(\frac{\dot{a}}{a})^2 = \frac{8\pi G}{3}\rho\tag{3.8}$$
这就是第一弗里德曼方程。带入\(H(t)= \dot{a}/a\)和\(\rho_{cr} = 3H^2_0/8\pi G\)，可以得到：
$$\frac{H^2(t)}{H^2_0} = \frac{\rho(t)}{\rho_{cr}}= \sum_{s=r,m,\nu,DE} \Omega_s [a(t)]^{-3(1+\omega_s)}\tag{3.9}$$
这里的能量密度包含了宇宙中所有物质，包括辐射、物质、中微子，暗能量。在我们上述的推导中，我们假设了宇宙是平坦的，所以式\((3.9)\)中还应该包含有宇宙曲率的部分。这部分非常简单，定义：\(\Omega_K \equiv 1-\Omega_0 \equiv 1- \sum \Omega_s\)，因此有：
$$\frac{H^2(t)}{H^2_0} = \sum_{s=r,m,\nu,DE} \Omega_s [a(t)]^{-3(1 + \omega_s)} + \Omega_K[a(t)]^{-2}\tag{3.10}$$
这就是均匀宇宙中物质演化的方程。
## 3.2 波尔兹曼方程
在相空间\(\{\bm{x},\bm{p}\}\)中，我们可以利用分布函数来描述粒子的数量：
$$N(\bm{x},\bm{p},t) = f(\bm{x},\bm{p},t) (\Delta x)^3 \frac{(\Delta p)^3}{(2\pi)^3} \tag{3.11}$$
现在我们想推导一个控制这个分布函数的方程。该方程应从单个粒子所服从的运动方程中唯一得出。我们忽略粒子之间的相互作用，则对粒子来说只有一个长程力作用，我们将其视作一个场：\(\bm{a}(\bm{x},\bm{p},t)\)。利用动量的定义，我们可以得到粒子的运动方程：
$$\bm{\dot{x}} = \frac{\bm{p}}{m};\bm{\dot{p}} = m \bm{a}(\bm{x},\bm{p},t) \tag{3.12}$$
由于粒子数守恒，所以分布函数的总时间导数为0：
$$\frac{d f(\bm{x},\bm{p},t)}{d t} = 0\tag{3.13}$$
带入运动方程，式\((3.13)\)变为：
$$\frac{\partial f}{\partial t} = -\frac{\bm{p}}{m} \nabla_x f - m\bm{a} \nabla_p f\tag{3.14}$$
这代表着一定相空间内，粒子数是守恒的；改变的粒子数是进入/离开相空间的粒子数。而当我们考虑粒子之间的相互碰撞时，需要在方程右侧添加一个碰撞项，用来描述粒子从一个相空间部分转移到另一个部分：
$$\frac{d f}{d t} = C[f]\tag{3.15}$$
完成上述推导后，我们可以从分布函数推导出所有宏观属性，如密度和压强等。比如给定分布函数时，能动张量的相对论表达式为：
$$T^\mu_\nu (\bm{x},t) = \frac{g}{\sqrt{-det[g_{\alpha\beta}]}} \int \frac{dP_1 dP_2 dP_3}{(2\pi)^3}\frac{P^\mu P_\nu}{P^0} f(\bm{x},\bm{p},t) \tag{3.16}$$
其中\(P^\mu\)是共动动量，\(\bm{p}\)为物理动量，\(g\)为简并因子。能动张量实际上给出了分布函数为\(f\)的粒子的动量流密度。
### 3.2.1 谐振子势下粒子的波尔兹曼方程
一维谐振子势下粒子的能量为：
$$E = \frac{p^2}{2m} + \frac{1}{2} k x^2\tag{3.17}$$
此时分布函数\(f\)为三个标量的函数：\(f(x,p,t)\)。假设粒子间没有相互作用，即式\(3.15\)中\(C[f]=0\)，则波尔兹曼方程为：
$$\frac{\partial f}{\partial t} + \frac{p}{m} \frac{\partial f}{\partial x} - kx \frac{\partial f}{\partial p} = 0\tag{3.18}$$
这里的第二项控制粒子在真实空间中移动的速度;前面的系数只是速度，\(v=p/m\)。最后一项控制粒子失去或获得动量的速度。

为了求解该方程，我们需要给定初始条件。但即使没有初始条件，我们也可以得到一个重要的结论：考虑一个静态分布，定义为\(\frac{\partial f}{\partial t} = 0\)（与 \(df/dt = 0 \)相反，后者始终成立）。这种分布，也称为平衡分布，意味着在空间\(x\)中的每个点，具有给定动量\(p\)的粒子数量在统计意义上保持不变。当然，这并不意味着粒子本身不会移动。对于该分布，一般解为：
$$f(p,x) = f_{EQ}(E[p,x]) = f_{EQ}(\frac{p^2}{2m} + \frac{1}{2} k x^2)\tag{3.19}$$
其中\(f_{EQ}\)是平衡分布函数。我们可以通过下面该式进行验证：
$$\frac{p}{m} \frac{\partial f(E)}{\partial x} - kx \frac{\partial f(E)}{\partial p} =\frac{d f}{d E}(\frac{p}{m} \frac{\partial E}{\partial x} - kx \frac{\partial E}{\partial p}) =0\tag{3.20}$$
### 3.2.2 膨胀宇宙中的波尔兹曼方程
在四维时空中，式\(3.12\)变为测地线方程，而三维空间动量变为四维动量：\(\P^\mu \equiv \frac{d x^\mu}{d \lambda}\)。然而，给定粒子集合的分布函数仍然是在六维相空间上定义的函数：首先，我们将时间分离出来；然后，每个粒子的四动量服从质壳约束：
$$P^2 \equiv g_{\mu\nu}P^\mu P^\nu = m^2\tag{3.21}$$
其中，\(P^2\)是四动量的平方，\(m\)是粒子的静质量。一般的三维动量为：\(p^2 \equiv g_{ij}P^iP^j\)，而四动量的时间分量为\(P^0 = \sqrt{p^2 + m^2}\)。因此，我们可以用\(p\)来表示\(P^0\)，所以分布函数还是\(f(\bm{x},\bm{p},t)\)，其中\(\bm{p}\)是三维动量。我们认为单位向量\(\hat{p}^i = p^i/p\)与\(P^i\)成正比：\(P^i \equiv C \hat{p}^i\)，其中\(C\)是一个与动量大小相关的系数。为了确定这个系数：
$$p^2 = g_{ij}\hat{p}^i\hat{p}^jC^2=a^2 \delta_{ij}\hat{p}^i\hat{p}^jC^2=a^2 C^2\tag{3.22}$$
由此我们得到\(C = p/a\)，所以：
$$P^i = \frac{p}{a} \hat{p}^i\tag{3.23}$$
所以我们重写式\(3.13\)：
$$\frac{d f}{d t} = \frac{\partial f}{\partial t} + \frac{\partial f}{\partial x^i} \frac{d x^i}{d t} + \frac{\partial f}{\partial p} \frac{d p}{d t} + \frac{\partial f}{\partial \hat{p}^i} \frac{d \hat{p}^i}{d t} = 0\tag{3.24}$$
与之前一样，我们在光滑膨胀宇宙中进行推导。正如我们在第二章中所见，粒子动量的方向\(\hat{p}^i\)在均匀宇宙中不会改变，所以最后一项为0。对于第二项，我们有：
$$\frac{d x^i}{d t} = \frac{P^i}{P^0} = \frac{p}{E} \frac{\hat{p}^i}{a}\tag{3.25}$$
对于第三项，我们需要计算\(\frac{d p}{d t}\)。我们从测地线方程开始：
$$\frac{d P^\mu}{d \lambda} + \Gamma^\mu_{\nu\rho} P^\nu P^\rho = 0\tag{3.26}$$
其时间分量为：
$$\frac{d P^0}{d \lambda} = -\Gamma^0_{\alpha\beta} P^\alpha P^\beta \tag{3.27}$$
对于\(dt/d\lambda = P^0\)的情况，我们有：
$$P^0\frac{d P^0}{d t} = -\Gamma^0_{ij} P^i P^j \tag{3.28}$$
由质壳约束，我们有：\(E^2 = (P^0)^2 = p^2 + m^2\)，所以：\(P^0 dP^0/d t = (1/2)d(E^2)/dt\)；而带入克氏符我们可以得到：
$$\frac{d p}{d t} = -H p\tag{3.29}$$
这个方程表明，在膨胀宇宙中，粒子的动量以\(1/a\)的速率衰减。我们可以将其带入式\(3.24\)中，得到：
$$\frac{\partial f}{\partial t} + \frac{p}{E} \frac{\hat{p}^i}{a} \frac{\partial f}{\partial x^i} - H p \frac{\partial f}{\partial p} = C[f]\tag{3.30}$$
这个式子对所有粒子都是适用的。然而，我们经常遇到两种极端情况：非相对论极限和相对论极限。
在相对论极限下，\(p\gg m\)，此时\(E\simeq p\)，由此我们得到：
$$\frac{\partial f}{\partial t} + \frac{1}{a} \hat{p}^i \frac{\partial f}{\partial x^i} - Hp \frac{\partial f}{\partial p} = C[f]\tag{3.31}$$
这适用于光子和中微子等相对论粒子。
在非相对论极限下，\(p\ll m\)，此时\(E\simeq m\)，所以：
$$\frac{\partial f}{\partial t} + \frac{p}{m a} \hat{p}^i \frac{\partial f}{\partial x^i} - H p \frac{\partial f}{\partial p} = C[f]\tag{3.32}$$
这适用于非相对论粒子，如电子和质子等。

在均匀宇宙中，\(\partial f/\partial x^i = 0\)，我们对式\(3.30\)进行积分，得到：
$$\int \frac{d^3 p}{(2\pi)^3}\frac{\partial f}{\partial t} - \int \frac{d^3 p}{(2\pi)^3}H p \frac{\partial f}{\partial p} =\int \frac{d^3 p}{(2\pi)^3} C[f]\tag{3.33}$$
左式第二项可以通过分部积分得到：
$$\int \frac{d^3 p}{(2\pi)^3} H p \frac{\partial f}{\partial p} = H \int \frac{d^3 p}{(2\pi)^3} (3f + p \frac{\partial f}{\partial p}) = 3H n\tag{3.34}$$
其中\(n = \int d^3 p/(2\pi)^3 f\)是粒子数密度。带入式\(3.33\)中，我们得到：
$$\dot{n} + 3H n = \int \frac{d^3 p}{(2\pi)^3} C[f]\tag{3.35}$$
这个方程描述了均匀宇宙中粒子数密度的演化。在没有碰撞的情况下，右侧为0，此时粒子数密度随时间的演化为：
$$n(t) = n_0 a^{-3}(t)\tag{3.36}$$
其中\(n_0\)是初始粒子数密度。这个结果表明，在均匀宇宙中，粒子数密度随着宇宙膨胀而减少，正比于\(a^{-3}\)。这个结果与我们在第二章中对非相对论粒子的讨论一致。
### 3.2.3 碰撞项
在玻尔兹曼领域，直接粒子相互作用的影响被称为“碰撞”。碰撞包括散射以及成对、湮灭和粒子衰减。一种典型的碰撞是两个例子相互作用后产生两个新粒子：
$$(1)_p + (2)_q \longleftrightarrow (3)_{p^\prime} + (4)_{p^\prime} \tag{3.37}$$
每种粒子都有独立的分布函数，\(f_s(\bm{x},\bm{p},t),s=1,2,3,4\)；但在宇宙学中，不同组分的分布函数通常是相同的，因此，我们不用单独的函数来跟踪它们，而是分配适当的统计权重\(g_s\)。
式\(3.37\)是如何影响粒子分布函数的演化的呢？首先，我们固定时间和空间，此时我们只需要讨论动量的影响。例如我们对粒子1进行讨论：我们需要减去失去的，并且加上获得的。失去的粒子数是\(f_1(\bm{p})\)乘以入射粒子1的动量分布函数\(f_2(\bm{q})\)，而获得的粒子数是\(f_3(\bm{p}^\prime)f_4(\bm{q}^\prime)\)。因此，碰撞项可以写为：
$$\begin{aligned}C[f_1(\bm{p})] =& \sum_{\bm{q},\bm{q}^\prime,\bm{p}^\prime}^{\bm{p}+\bm{q} = \bm{p}^\prime + \bm{q}^\prime} \delta^{(1)}_D (E_1(p) + E_2(q)- E_3(p^\prime)- E_4(q^\prime))\\&|\mathcal{M}|^2 \times \{f_3(\bm{p}^\prime)f_4(\bm{q}^\prime ) - f_1(\bm{p}) f_2(\bm{q})\}\tag{3.38}\end{aligned}$$
其中\(\mathcal{M}\)是散射振幅，\(\delta^{(1)}_D\)是狄拉克δ函数，表示能量守恒。在这里，散射幅度的平方\(|\mathcal{M}|^2\)取决于相互作用的微观物理细节，可以使用费曼图进行计算。这个方程表明，碰撞项的变化率与入射粒子和出射粒子的分布函数有关。
在上面的推导中，我们忽略了如泡利法则等量子效应；加入之后，会在正反应过程中加入\((1 \pm f_3)(1 \pm f_4)\)的因子，逆反应过程加入\((1 \pm f_1)(1 \pm f_2)\)的因子。当粒子为玻色子时，为加号；为费米子时，用减号。
由狄拉克函数的复合函数公式\(\delta^{(1)}_D(f(x)) = \sum_i \frac{\delta^{(1)}_D(x-x_i)}{|f'(x_i)|}\)，我们有：
$$\int d^3p\int_{0}^{\infty}dE\delta_{\mathrm{D}}^{(1)}(E^2-p^2-m^2)=\int d^3p\int_{0}^{\infty}dE\frac{\delta_{\mathrm{D}}^{(1)}\left(E-\sqrt{p^2+m^2}\right)}{2E}\tag{3.39}$$
因此，碰撞项可以写为：
$$\begin{aligned}
C[f_1(p)] & 
\begin{aligned}
 & =\frac{1}{2E_1(p)}\int\frac{d^3q}{(2\pi)^32E_2(q)}\int\frac{d^3p^{\prime}}{(2\pi)^32E_3(p^{\prime})}\int\frac{d^3q^{\prime}}{(2\pi)^32E_4(q^{\prime})}|\mathcal{M}|^2
\end{aligned} \\
 & \times(2\pi)^4\delta_\mathrm{D}^{(3)}[p+q-p^{\prime}-q^{\prime}]\delta_\mathrm{D}^{(1)}[E_1(p)+E_2(q)-E_3(p^{\prime})-E_4(q^{\prime})] \\
 & \times\left\{f_3(\boldsymbol{p}^{\prime})f_4(\boldsymbol{q}^{\prime})\Big[1\pm f_1(\boldsymbol{p})\Big]\Big[1\pm f_2(\boldsymbol{q})\Big]\right. \\
 & -f_1(\boldsymbol{p})f_2(\boldsymbol{q})\Big[1\pm f_3(\boldsymbol{p}^{\prime})\Big]\Big[1\pm f_4(\boldsymbol{q}^{\prime})\Big]\Big\}.
\end{aligned}\tag{3.40}$$
这个方程描述了粒子分布函数在碰撞过程中的变化。它表明，粒子分布函数的变化率与入射和出射粒子的分布函数以及散射振幅的平方有关。
## 3.3 非均匀宇宙
### 3.3.1 时空微扰
我们选取FLRW度规作为背景度规，引入微扰\(\Psi\)和\(\Phi\)：
$$\begin{aligned}
 & g_{00}(x,t)=-1-2\Psi(x,t), \\
 & g_{0i}(x,t)=0, \\
 & g_{ij}(x,t)=a^2(t)\delta_{ij}\left[1+2\Phi(x,t)\right].
\end{aligned}\tag{3.41}$$
在没有微扰的情况下，这是FLRW度规，描述的是均匀各向同性的宇宙；而若是\(a(t)=1\)，则是弱场极限下的牛顿引力势\(\Psi\)和空间曲率\(\Phi\)。对于(3.41)中的微扰，我们有两个技术点：
- 首先，可以将扰动分解为在从一个 3D 空间坐标系到另一个 3D 空间坐标系的变换下表现为标量、向量和张量的扰动。这个式子中的\(\Psi\)和\(\Phi\)是标量扰动，描述了时空的标量扰动；而在更高阶的微扰中，我们还可以引入向量和张量扰动。
- 其次，张量微扰可以与引力波联系起来。

值得注意的式(3.41)的另一个特征是它的形式对应于度规或规范的特定选择。与电磁相互作用一样，同样的场可以选取不同的规范。我们将式(3.41)称作共形牛顿规范。
由于\(\Psi\)和\(\Phi\)是微扰量，我们只考虑一价量级的项。我们可以计算克里斯托弗符：
$$\begin{aligned}
 & \Gamma^0_{00}=\dot{\Psi}, \\
 & \Gamma^0_{0i}=\Psi_{,i}, \\
 & \Gamma^0_{ij}=a^2\left[\delta_{ij}(H-2H\Phi+\dot{\Phi})-2\Phi_{,ij}\right], \\
 & \Gamma^i_{00}=\frac{1}{a^2}\Psi_{,i}, \\
 & \Gamma^i_{0j}=(H+\dot{\Phi})\delta_{ij}, \\
 & \Gamma^i_{jk}=\Phi_{,k}\delta_{ij}+\Phi_{,j}\delta_{ik}-\Phi_{,i}\delta_{jk}.
\end{aligned}\tag{3.42}$$
由此，我们可以计算里奇张量：
$$\begin{aligned}
 & R_{00}=3\dot{H}-3\ddot{\Phi}-3H\dot{\Psi}+3H\dot{\Phi}-\frac{1}{a^2}\nabla^2\Psi, \\
 & R_{0i}=2\Psi_{,i}, \\
 & R_{ij}=\delta_{ij}\left[-(2\dot{H}+3H^2)+2\ddot{\Phi}+6H\dot{\Phi}-2H\dot{\Psi}+\frac{2}{3a^2}\nabla^2\Phi\right] \\
 & -\frac{1}{a^2}\left(\Phi_{,ij}-\frac{1}{3}\delta_{ij}\nabla^2\Phi+\Psi_{,ij}-\frac{1}{3}\delta_{ij}\nabla^2\Psi\right).
\end{aligned}\tag{3.43}$$
由此，我们可以计算里奇标量：
$$R = 6(\dot{H}+2H^2) - 6\ddot{\Phi} - 12H\dot{\Phi} + 6H\dot{\Psi} + \frac{2}{a^2}\nabla^2\Psi -\frac{4}{a^2}\nabla^2\Phi.\tag{3.44}$$
### 3.3.2 测地线方程
为了推导波尔兹曼方程，我们需要知道粒子如何在微扰时空中运动，即我们需要得到测地线方程。这个过程中我们最终需要的是\(\frac{d \hat{x}^i}{d t}\)、\(\frac{d p}{d t}\)和\(\frac{d \hat{p}^i}{d t}\)。在微扰时空中的质壳约束为：
$$g_{\mu\nu}P^\mu P^\nu=-(1+2\Psi)(P^0)^2+p^2=-m^2\tag{3.45}$$
此时同样有\(p^2 \equiv g_{ij}P^i P^j\)和\(E(p) \equiv \sqrt{p^2 + m^2}\)。

对于无质量粒子，动量的时间分量为：
$$P^0=\frac{E}{\sqrt{1+2\Psi}}=E(1-\Psi).\tag{3.46}$$
同理，我们对空间分量进行计算，最终得到：
$$P^\mu=\left[E(1-\Psi),p^i\frac{1-\Phi}{a}\right].\tag{3.47}$$
与之前相同，我们同样有：
$$p^i=p\hat{p}^i\quad\mathrm{where}\quad\hat{p}^i=\hat{p}_i \quad and \quad\delta_{ij}\hat{p}^i \hat{p}^j =   0\tag{3.48}$$
之前我们定义了\(P^i \equiv d x^i/d\lambda\)和\(P^0 \equiv d t / d \lambda\)，所以我们有：
$$\begin{aligned}
\frac{dx^i}{dt} & =\frac{dx^\iota}{d\lambda}\frac{d\lambda}{dt} \\
 & =\frac{P^i}{P^0}=\frac{\hat{p}^i}{a}\frac{p}{E}\left(1-\Phi+\Psi\right).
\end{aligned}\tag{3.49}$$
为了计算\(\frac{d p}{d t}\)和\(\frac{d \hat{p}^i}{d t}\)我们需要先得到\(\frac{d \hat{p}^i}{d t}\)。我们可以利用测地线方程：
$$\begin{aligned}
\frac{dp^i}{d\lambda} & =\frac{d}{d\lambda}\left[(1+\Phi)aP^i\right] \\
 & =P^i\frac{d}{d\lambda}[(1+\Phi)a]+(1+\Phi)a\frac{dP^i}{d\lambda}.
\end{aligned}\tag{3.50}$$
第一项可以利用\(d/d\lambda = P^\mu \partial/\partial x^\mu\)，得到：
$$\frac{d}{d\lambda}[(1+\Phi)a]=P^0a[H+\dot{\Phi}]+aP^k\Phi_{,k}.\tag{3.51}$$
第二项可以利用测地线方程得到：
$$\begin{aligned}
\frac{dP^i}{d\lambda} & =-\Gamma^i_{\alpha\beta} P^\alpha P^\beta \\
 & =-\left[\Gamma^i{}_{00}P^0P^0+2\Gamma^i{}_{0j}P^0P^j+\Gamma^i{}_{jk}P^jP^k\right].
\end{aligned}\tag{3.52}$$
代入式（3.42）中的克氏符，有：
$$\frac{dP^i}{d\lambda}=-E\left\{\frac{E}{a^2}\Psi_{,i}+2\left(H+\dot{\Phi}\right)\frac{p^i}{a}(1-\Psi-\Phi)+\frac{2}{a^2}\frac{p^i}{E}p^k\Phi_{,k}-\frac{p^2}{a^2E}\Phi_{,i}\right\}.\tag{3.53}$$
因此，我们有：
$$\begin{aligned}
\frac{dp^i}{d\lambda} & =E(1-\Psi)\left\{[H+\dot{\Phi}]p^i+p^k\Phi_{,k}\frac{p^i}{aE}\right\} \\
 & -E\left\{\frac{E}{a}\Psi_{,i}+2\left(H+\dot{\Phi}\right)p^i(1-\Psi)+\frac{2}{a}\frac{p^i}{E}p^k\Phi_{,k}-\frac{p^2}{aE}\Phi_{,i}\right\}.
\end{aligned}\tag{3.54}$$
代入\(dp^i/dt=(P^0)^{-1}dp^i/d\lambda\)：
$$\frac{dp^i}{dt}=-\left(H+\dot{\Phi}\right)p^i-\frac{E}{a}\Psi_{,i}-\frac{1}{a}\frac{p^i}{E}p^k\Phi_{,k}+\frac{p^2}{aE}\Phi_{,i}.\tag{3.55}$$
代入\(\frac{dp}{dt}=\frac{d}{dt}\sqrt{\delta_{ij}p^ip^j}=\delta_{ij}\frac{p^i}{p}\frac{dp^j}{dt}\)：
$$\begin{aligned}
\frac{dp}{dt} & =-\left[H+\dot{\Phi}\right]p-\frac{E}{a}\hat{p}^{i}\Psi_{,i}-\frac{1}{a}\frac{p^{2}}{E}\hat{p}^{k}\Phi_{,k}+\frac{p^{2}}{aE}\hat{p}^{i}\Phi_{,i} \\
 & =-\left[H+\dot{\Phi}\right]p-\frac{E}{a}\hat{p}^i\Psi_{,i}.
\end{aligned}\tag{3.56}$$
这个式子描述了在FLRW度规下的微扰宇宙中，粒子动量的演化。它表明，粒子的动量随着宇宙膨胀而衰减，并且受到时空微扰的影响。第一项解释了由于哈勃膨胀而失去的动量，我们将\(H+\dot{\Phi}\)定义为本地膨胀率；第二项描述了粒子在引力场下的动量变化：进入引力场时，粒子动量增加；离开引力场时，粒子动量减少。这与牛顿引力场下的粒子动量变化一致。
注意到式（3.55）中，\(\dot{\Phi}_{,i}\)和\(p^k  \Phi_{,k}\)并没有出现在（3.56）中。它们不会改变粒子的线性顺序动量，但它们确实会改变动量的方向。为了证明这一点，我们可以计算得到：
$$\begin{aligned}
\frac{d\hat{p}^i}{dt} & =\frac{1}{p}\frac{dp^i}{dt}-\frac{p^i}{p^2}\frac{dp}{dt} \\
 & =\frac{E}{ap}\left[\delta^{ik}-\hat{p}^i\hat{p}^k\right]\left(\frac{p^2}{E^2}\Phi-\Psi\right)_{,k}.
\end{aligned}\tag{3.57}$$
势中的空间梯度改变了大质量粒子和无质量粒子的轨迹，而测地线方程代表弯曲空间中的直线，这说明引力势场改变了时空的结构。
### 3.3.3 辐射的无碰撞波尔兹曼方程
将（3.49）和（3.55）代入波尔兹曼方程（3.30），我们得到：
$$\begin{aligned}
\frac{df}{dt} & =\frac{\partial f}{\partial t}+\frac{\partial f}{\partial x^i}\frac{\hat{p}^i}{a}\left(1-\Phi+\Psi\right)-\frac{\partial f}{\partial p}\left\{\left[H+\dot{\Phi}\right]p+\frac{1}{a}p^i\Psi_{,i}\right\} \\
 & +\frac{\partial f}{\partial\hat{p}^i}\frac{1}{a}\left[(\Phi-\Psi)_{,i}-\hat{p}^i\hat{p}^k\left(\Phi-\Psi\right)_{,k}\right].
\end{aligned}\tag{3.58}$$
当然，我们同样忽略高阶项，即只保留一价项。此时我们需要注意到，\(\frac{\partial f}{\partial\hat{p}^i}\)和\(\frac{\partial f}{\partial x^i}\)都属于一价项（因为零价项代表的是均匀宇宙，与方向无关），所以我们可以将式（3.58）写成：
$$\frac{df}{dt}=\frac{\partial f}{\partial t}+\frac{\hat{p}^i}{a}\frac{\partial f}{\partial x^i}-\left[H+\dot{\Phi}+\frac{1}{a}\hat{p}^i\frac{\partial\Psi}{\partial x^i}\right]p\frac{\partial f}{\partial p}.\tag{3.59}$$
式 （3.59） 将直接引导我们了解控制 CMB 各向异性的方程。
### 3.3.4 带质量粒子的无碰撞波尔兹曼方程
与辐射类似，我们将（3.49）和（3.56）代入波尔兹曼方程（3.32），得到：
$$\begin{gathered}
\frac{df}{dt}=\frac{\partial f}{\partial t}+\frac{\partial f}{\partial x^i}\frac{\hat{p}^i}{a}\frac{p}{E}\left(1-\Phi+\Psi\right)-p\frac{\partial f}{\partial p}\left[H+\dot{\Phi}+\frac{E}{ap}\hat{p}^i\Psi_{,i}\right] \\
+\frac{\partial f}{\partial\hat{p}^i}\frac{E}{ap}\left[\left(\frac{p^2}{E^2}\Phi-\Psi\right)_{,i}-\hat{p}^i\hat{p}^k\left(\frac{p^2}{E^2}\Phi-\Psi\right)_{,k}\right].
\end{gathered}\tag{3.60}$$
同样，我们忽略高阶项，得到：
$$\frac{df}{dt}=\frac{\partial f}{\partial t}+\frac{p}{E}\frac{\hat{p}^i}{a}\frac{\partial f}{\partial x^i}-\left[H+\dot{\Phi}+\frac{E}{ap}\hat{p}^i\Psi_{,i}\right]p\frac{\partial f}{\partial p}.$$
这与辐射的情况类似，但有质量粒子的影响。