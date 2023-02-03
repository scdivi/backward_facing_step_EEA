# backward_facing_step_EEA

This script compares the rate of convergence for a quantity of interest with both [uniform](https://github.com/scdivi/backward_facing_step_EEA/tree/uniform)
and [residual-based](https://github.com/scdivi/backward_facing_step_EEA) adaptive refinements using a known problem called backward facing step.

## Problem definition

![image](https://user-images.githubusercontent.com/33148729/216629499-95361905-7bcb-4420-8ad7-5aa334b8f2a8.png)

### Strong form
$$
			\begin{equation*}
				\begin{aligned}
					&\text{Given constants}\, \mu \, \text{and}\, \gamma, &&\\
					&\text{find}\, \mathbf{u}: \Omega \in \mathbb{R}^2\, \text{such that:} && \\
					& - {\rm div}(\boldsymbol{\sigma}) = 0 && \text{in}\, \Omega\\
					&\phantom{-div(u)}\mathbf{u} = \mathbf{0}   && \text{on}\, \Gamma_{D} \\
					&\phantom{-div(u} u_{x} = u_{\rm max}(x_1 - h) \cfrac{H - (x_1 - h)}{\left(\frac{H}{2}\right)^2} && \text{on}\, \Gamma_{in} \\
					&\phantom{-div} \boldsymbol{\sigma} \cdot \mathbf{n} = \mathbf{0} && \text{on}\, \Gamma_{out} \\
				\end{aligned}
			\end{equation*}
$$

with

$$
		\begin{align*}
			\boldsymbol{\sigma} = \mu (2 \boldsymbol{\varepsilon} + \gamma {\rm tr}(\boldsymbol{\varepsilon}) \mathbf{I}) \quad \quad
			\boldsymbol{\varepsilon} = \frac{1}{2} \left( \nabla \mathbf{u} + (\nabla \mathbf{u})^{\rm T} \right)
		\end{align*}
$$

### Weak form

$$
			\begin{equation*}
				\begin{aligned}
					&\text{Find}\, \mathbf{u} \in \mathcal{S} \text{such that:} &&\\
					& \mathcal{B}(\mathbf{u},\mathbf{v}) = \mathcal{F} (\mathbf{v}) && \forall \mathbf{v} \in \mathcal{V}
				\end{aligned}
			\end{equation*}
$$

with

$$
			\begin{align*}
				\mathcal{B}(\mathbf{u},\mathbf{v}) &= \int_{\Omega} \left[ 2 \mu \boldsymbol{\varepsilon} \left( \frac{1}{2} \left( \nabla \mathbf{v} + (\nabla \mathbf{v})^{\rm T} \right) \right)  + \gamma \mu \, {\rm div}(\mathbf{u}) \, {\rm div}(\mathbf{v}) \right]\, {\rm dV}\\
				\mathcal{F}(\mathbf{v}) &= \int_{\Gamma_{out}} \mathbf{v} \cdot (\boldsymbol{\sigma} \cdot \mathbf{n}) \, {\rm dS} = 0
			\end{align*}
$$

### Residual-based error

$$
\begin{equation*}
					\mathcal{R}^h(\mathbf{v}) := B(\mathbf{u}^h, \mathbf{v}) - \mathcal{F}(\mathbf{v}) \phantom{ {\sup_{\mathbf{v} \in \mathcal{V}\setminus 0}} \frac{\mathcal{R}^H(\mathbf{v})}{\| \mathbf{v} \|_{\mathcal{V}}}}
				\end{equation*}	
$$
