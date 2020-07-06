# 数字信号处理

## 绪论

### 基本概念

#### 信号：

* 概念：信息的物理表现方式/传递信息的函数
	* 一维信号：函数中只有一个独立变量
		* 语音信号、温度信号等
	* 二维信号：函数中有两个独立变量表示
		* 一副黑白图像：空间坐标
	* 多维信号：函数优多个独立变量表示
		* 一副彩色图像：空间座标
* 分类
	* 周期信号/非周期信号
	* 确定信号/随机信号
	* 连续信号（模拟信号）/时域离散信号/幅度离散信号/数字信号

#### 系统

* 概念：将信号进行处理（或变换），以达到人们要求的各种设备。
* 分类
	* 连续信号（模拟信号）系统
	* 时域离散信号系统
	* 幅度离散信号系统
	* 数字信号系统

#### 信号处理

* 概念：研究用系统对含有信息的信号进行处理（变换），获得人们所希望的信号，提取信息，标语利用的一门学科

* 内容

	滤波	变换	分析	估值	检测	压缩	识别

* 分类

	* 模拟信号处理：模拟期间
	* 数字信号处理（DSP）：数值计算

#### 数字信号处理

* 概念：利用寒粟子计算机或专用数字（集成）电路按预订规则对数字信号进行数值运算

#### DSP的实现方法

* 软件实现的方法
	* 缺点：运算速度慢，达不到实时
* 硬件实现方法
	* 缺点：不灵活
* 软硬结合法
	* 特点：既灵活又快速，工程技术领域中的主要实现方法

## 第一部分 数字信号的表示和分析工具

### 第一章

#### 1.1 时域离散信号

##### 序列

1. 定义：对模拟信号$x_a(t)$进行等间隔采样，采样间隔为 T ，得到离散时间信号（序列）
	$$
	x(n)=x_a(t) \mid_{t=nT}=x_a(nT), \quad \text{$ -\infty < n<+ \infty$}
	$$

2. 表示方式

  * 集合表示
$$
x(n)=\{ \dots 1,3,2.5,3.3,1.9,0,4.1 \dots \}
$$

  * 公式表示
$$
  	x(n)=x_a(nT), \quad \text{$ -\infty < n < + \infty $} \\
  	x(n)=x_a(t) \mid_{t=nT}=0.9 sin(50 \pi nT)
$$

  * 图形表示

3. 常用的典型序列

  * 单位采样序列
$$
  	\delta(n)=
  	\begin{cases}
  	1, & \text{$ n=0 $} \\
  	0, & \text{$ n \not= 0 $}
  	\end{cases}
$$
  	特点：尽在$n=0$时取值为1，其他均为零。

  * 单位阶跃序列
$$
  	u(n)=
  	\begin{cases}
  	1, & \text{$n \geq 0$} \\
  	0, & \text{$n < 0$}
  	\end{cases}
$$
  	与单位采样序列的关系：
$$
  	\delta(n)=u(n)-u(n-1) \\
  	
  	\begin{aligned}
  	u(n) & =\sum_{k=0}^{\infty} \delta(n) + \delta(n-1) + \delta(n-2) + \cdots \\
  & =\sum_{m=-\infty}^n \delta(m) , \quad \text{令$n-k=m$}
  	\end{aligned}
$$

  * 矩形序列
$$
  	R_N(n)=
  	\begin{cases}
  	1 & \text{$0 \leq n \leq N-1$} \\
  	0 &	\text{其他$n$}
  	\end{cases}
$$
  	与其他序列的关系
$$
  	R_N(n)=u(n)-u(n-N) \\
  	R_N(n)=\sum_{m=0}^{N-1}\delta(n-m)=\delta(n)+\delta(n-1)+\cdots+\delta[n-(N-1)]
$$


  * 实指数序列
$$
  	\begin{aligned}
  	& x(n)=a^n u(n) \quad \text{$a$ 为实数} \\
  	& |a|<1,x(n)的幅度随n的增大而减小（收敛） \\
  	& |a|>1,发散 
  	\end{aligned}
$$


  * 正弦序列
$$
  	x(n)=sin(\omega n)
$$
  	若正弦序列是由模拟正弦信号采样得到，则
$$
  	\begin{aligned}
  	& x_a(t)=sin(\Omega t) \\
  	& x_a(t)\mid_{t=nT}=sin(\Omega nT) \\
  	& x(n)=sin(\omega n)
  	\end{aligned}
$$
  	数字频率 $\omega$ 与模拟角频率 $\Omega$ 关系：
$$
  	\omega = \Omega T \\
  	\omega = \frac{\Omega}{F_s}
$$
  	其中：
$$
  	\begin{aligned}
  	& \omega: 数字角频率 \\
  	& \Omega: 模拟角频率 \\
  	& T: 采样间隔/周期 \\
  	& F_s: 采样频率
  	\end{aligned}
$$


  * 复指数序列
$$
  	\begin{aligned}
  	x(n) & = e^{(\sigma+j \omega_0)n} = e^{\sigma n} \cdot e^{j \omega_0 n} \\
  	& = e^{\sigma n} cos(\omega_0 n)+je^{\sigma n} sin(\omega_0 n) \\
  	& -其中 \omega_0 为数字域角频率
  	\end{aligned}
$$

$$
  	\begin{aligned}
  	
  	& e^{j \omega_0 n}=cos(\omega_0 n)+ j sin(\omega_0 n) \\
  	e^{j(\omega_0 + 2 \pi M)n} & = cos[(\omega_0 2\pi M)n]+jsin[(\omega_0+2 \pi M)n] \\ 
  	& =cos(\omega_0 n) + j sin(\omega_0 n) = e^{j \omega_0 n}
  	
  	\end{aligned}
$$

  	欧拉公式：
$$
\begin{aligned}
  	e^{j \omega} & = cos \omega + j sin \omega \\
  	\implies & e^{j \omega} + e^{-j \omega}=2cos \omega \\
  	& e^{j \omega} - e^{-j \omega} = 2j sin \omega
  	\end{aligned}
$$
  	正弦序列和复指数序列均以 $2\pi$ 为周期，频率域分析时只研究主值区间即可。 

  * 周期序列

  	如果对于所有 $n$ 存在一个最小的正整数 N ，满足：
$$
  	x(n)=x(n+N) \quad \quad  \text{$-\infty < n < \infty$}
$$
  	则称序列 $x(n)$ 是周期序列，周期为 N。

  	例：
$$
  	x(n)=sin(\frac{\pi}{4} n) \\
  	\omega = \frac{\pi}{4} \\
$$
  	由于 n 取整数，则:
$$
x(n)=sin(\frac{\pi}{4}(n+8)) \quad \quad x(n)是周期为8的周期序列
$$


  * 有限长任意序列

  	$x(n)$ 可以表示成单位采样序列的移位加权和，也可以表示成与单位采样序列的卷积和。
$$
x(n)=\sum_{m= -\infty}^{\infty}x(m) \delta(n-m)=x(n) * \delta(n)
$$

  

  













## 第二部分 数字滤波器的设计

 