# The origin of species

The very early universe was hot and dense, so the interactions of particles ocurred more frequently than in the present-day universe. These multiple interactions kept many of the constituents in the universe in equilibrium. Nonetheless, there were times when reactions could not proceed rapidly enough to maintain equilibrium conditions.

Indeed, we will see in this chapter that out-of-equilibrium phenomena played a role in

- (i) the formation of the light elements during Big Bang Nucleosynthesis;
- (ii) recombination of electrons and protons into neutral hydrogen; and possibly in
- (iii) the production of dark matter in the early universe.

## 4.1 The homogeneous Boltzmann equation revisited

Let us recall the interaction of two particles:\(1+2 \rightleftharpoons 3+4\), and when we are interested in the number density \(n_1\) of particle 1, we could use this two equations:

$$
\frac{dn(t)}{dt}+3Hn(t)=\int\frac{d^3p}{(2\pi)^3}C[f],\tag{1}
$$

$$
\begin{aligned}
C[f_1(p)] & =\frac{1}{2E_1(p)}\int\frac{d^3q}{(2\pi)^32E_2(q)}\int\frac{d^3p^{\prime}}{(2\pi)^32E_3(p^{\prime})}\int\frac{d^3q^{\prime}}{(2\pi)^32E_4(q^{\prime})}|\mathcal{M}|^2 \\
 & \times(2\pi)^4\delta_\mathrm{D}^{(3)}[p+q-p^{\prime}-q^{\prime}]\delta_\mathrm{D}^{(1)}[E_1(p)+E_2(q)-E_3(p^{\prime})-E_4(q^{\prime})] \\
 & \times\Big\{f_3(p^{\prime})f_4(q^{\prime})\Big[1\pm f_1(p)\Big]\Big[1\pm f_2(q)\Big] \\
 & -f_1(\boldsymbol{p})f_2(\boldsymbol{q})\Big[1\pm f_3(\boldsymbol{p}^{\prime})\Big]\Big[1\pm f_4(\boldsymbol{q}^{\prime})\Big]\Big\}.
\end{aligned}\tag{2}
$$

Combining Eq.(1) and Eq.(2), we obtain:

$$
\begin{aligned}
a^{-3}\frac{d(n_1a^3)}{dt} & =\int\frac{d^3p_1}{(2\pi)^32E_1}\int\frac{d^3p_2}{(2\pi)^32E_2}\int\frac{d^3p_3}{(2\pi)^32E_3}\int\frac{d^3p_4}{(2\pi)^32E_4} \\
 & \times(2\pi)^4\delta_\mathrm{D}^{(3)}(\boldsymbol{p}_1+\boldsymbol{p}_2-\boldsymbol{p}_3-\boldsymbol{p}_4)\delta_\mathrm{D}^{(1)}(E_1+E_2-E_3-E_4)|\mathcal{M}|^2 \\
 & \times\{f_3f_4[1\pm f_1][1\pm f_2]-f_1f_2[1\pm f_3][1\pm f_4]\}.
\end{aligned}\tag{3}
$$

We have thus obtained an integrodifferential equation for the phase-space distributions. Further, in principle at least, it must be supplemented with similar equations for the other species. However, in practice, these difficulties could be overcome  for many cosmology applications.
Scattering processes tipycally enforce kinetic equilibrium, leading to that the distributions of the various species take on the generic Bose–Einstein/Fermi–Dirac forms with equal temperature \(T\) for each species. This form condenses all of the freedom in the distribution into the functions of time \(T\) and \(\mu\). In the out-of-equilibrium cases we will study, the system will not be in chemical equilibrium and we will have to solve a differential equation for \(\mu\). And we could just solve the ODE not the complicated integral equation.

Considering the situation where \(E-\mu \gg 1\), Bose–Einstein and Fermi–Dirac distributions can approximately become Boltzmann distribution(which means ignoring the complications of quantum statistics):

$$
f(E)\to e^{\mu/T}e^{-E/T}.\tag{4}
$$

In addition, the Pauli blocking/Bose enhancement factors in the Boltzmann equation can be neglected. This is because \(f(E) \ll 1\) in the limit of Eq.(4).
Under these approximations, the last line of Eq. (3) becomes:

$$
\begin{aligned}
f_3f_4[1\pm f_1][1\pm f_2] & -f_1f_2[1\pm f_3][1\pm f_4] \\
 & \to\quad e^{-(E_1+E_2)/T}\left\{e^{(\mu_3+\mu_4)/T}-e^{(\mu_1+\mu_2)/T}\right\}.
\end{aligned}\tag{5}
$$

Here we have used energy conservation,\(E_1+E_2=E_3+E_4\), and the fact that the number of particles is conserved.We will use the number densities themselves as the time-dependent functions to be solved for, instead of\(\mu\). The number density of species \(s\) is related to \(\mu_s\) via

$$
n_s=g_se^{\mu_s/T}\int\frac{d^3p}{(2\pi)^3}e^{-E_s(p)/T},\tag{6}
$$

if \(\mu_s = 0\), then the equation:

$$
\left.n_s^{(0)}\equiv g_s\int\frac{d^3p}{(2\pi)^3}e^{-E_s(p)/T}=\left\{
\begin{array}
{ll}g_s\left(\frac{m_sT}{2\pi}\right)^{3/2}e^{-m_s/T} & m_s\gg T, \\
g_s\frac{T^3}{\pi^2} & m_s\ll T.
\end{array}\right.\right.\tag{7}
$$

With this definition, \(e^{\mu_i/T}\)can be rewritten as \(n_i/n_i^{(0)}\), and the last line of Eq. (5) becomes:

$$
e^{(\mu_3+\mu_4)/T}-e^{(\mu_1+\mu_2)/T}=\frac{n_3n_4}{n_3^{(0)}n_4^{(0)}}-\frac{n_1n_2}{n_1^{(0)}n_2^{(0)}}.\tag{8}
$$

With these approximations, the Boltzmann equation now simplifies enormously. Define the thermally averaged cross section as

$$
\begin{gathered}
\left\langle\sigma\upsilon\right\rangle
\begin{aligned}
 & \equiv\frac{1}{n_1^{(0)}n_2^{(0)}}\int\frac{d^3p_1}{(2\pi)^32E_1}\int\frac{d^3p_2}{(2\pi)^32E_2}\int\frac{d^3p_3}{(2\pi)^32E_3}\int\frac{d^3p_4}{(2\pi)^32E_4}e^{-(E_1+E_2)/T}
\end{aligned} \\
\times(2\pi)^4\delta_\mathrm{D}^{(3)}(\boldsymbol{p}_1+\boldsymbol{p}_2-\boldsymbol{p}_3-\boldsymbol{p}_4)\delta_\mathrm{D}^{(1)}(E_1+E_2-E_3-E_4)|\mathcal{M}|^2,
\end{gathered}\tag{9}
$$

which in general depends on the temperature T . Then, the Boltzmann equation becomes

$$
a^{-3}\frac{d(n_1a^3)}{dt}=n_1^{(0)}n_2^{(0)}\langle\sigma v\rangle\left\{\frac{n_3n_4}{n_3^{(0)}n_4^{(0)}}-\frac{n_1n_2}{n_1^{(0)}n_2^{(0)}}\right\}.\tag{10}
$$

One qualitative note about Eq. (10). The left-hand side is of order \(n_1/t\), or, since the typical cosmological time is \(H^{-1}\), \(n_1 H\) . The right-hand side is of order \(n_1n_2\langle\sigma v\rangle\). Therefore, if the reaction rate for a single particle of type 1,\(n_2\langle\sigma v\rangle\), is much larger than the expansion rate, then the terms on the right side will be much larger than the one on the left. The only way to maintain equality then is for the individual terms on the right to cancel. Thus, when reaction rates are large, Eq. (10) approaches

$$
\frac{n_3n_4}{n_3^{(0)}n_4^{(0)}}=\frac{n_1n_2}{n_1^{(0)}n_2^{(0)}}.\tag{11}
$$

This equation is equivalent to the condition \(\mu_1+\mu_2=\mu_3+\mu_4\), the relation we have referred to as chemical equilibrium above. The same relation is also known under the names of nuclear statistical equilibrium and Saha equation.

## 4.2 Big Bang nucleosynthesis

Big Bang Nucleosynthesis (BBN) happened when the temperature of the universe cooled to 1 MeV. At that point in time, the cosmic plasma consisted of:

- Relativistic particles in equilibrium: photons, electrons and positrons.Besides a small difference due to fermion/boson statistics, these all had the same abundances.
- Decoupled relativistic particles: neutrinos.Neutrinos share the same temperature as the other relativistic particles (but see Sect. 2.4.4), and hence are roughly as abundant, but they do not couple to them.
- Nonrelativistic particles: baryons.If there had been no asymmetry in the initial number of baryons and anti-baryons, then both would be completely depleted by 1 MeV. However, the asymmetry has happened. Comparing the abundance of baryons to photons, we find \(n_b/s \sim 10^{-10}\) today. There are many fewer baryons than relativistic particles in the universe.

To solve the Eq.(10) and at least for a qualitative understanding of the result, we can make use of two simplifications that obviate the need to solve the full set of differential equations.

- no elements heavier than helium are produced at appreciable levels. So the only nuclei that we need to trace are hydrogen and helium, and their isotopes: deuterium, tritium, and \(^3 He\).
- even in the context of this reduced set of elements, the physics splits up neatly into two parts since above \(T \simeq 0.1 MeV\), no light nuclei form: only free protons and neutrons exist.

Therefore, we first solve for the neutron/proton ratio and then use this abundance as input for the synthesis of helium and isotopes such as deuterium.
However, both of these simplifications rely on the physical fact that, at high temperatures, comparable to nuclear binding energies, any time a nucleus is produced in a reaction, it is destroyed by a high-energy photon.
This fact is reflected in the fundamental equilibrium equation (11).Let us consider the reaction:\(n + p \leftrightarrow D+ \gamma\). Since photons have \(n_\gamma = n_\gamma^{(0)}\), we can write:

$$
\frac{n_\mathrm{D}}{n_nn_p}=\frac{n_\mathrm{D}^{(0)}}{n_n^{(0)}n_p^{(0)}}.\tag{12}
$$

Using Eq. (7) for the quantities on the right-hand side leads to

$$
\frac{n_\mathrm{D}}{n_nn_p}=\frac{3}{4}\left(\frac{2\pi m_\mathrm{D}}{m_nm_pT}\right)^{3/2}e^{[m_n+m_p-m_\mathrm{D}]/T},\tag{13}
$$

the factor of 3/4 being due to the number of spin states (3 for D and 2 each for p and n). In the prefactor, \(m_D\) could be set to \(2m_n = 2 m_p\); and in the exponent, we could use the fact that \(m_n + m_p - m_D = B_D = 2.22MeV\). Therefore, we have:

$$
\frac{n_{\mathrm{D}}}{n_{n}n_{p}}=\frac{3}{4}\left(\frac{4\pi}{m_{p}T}\right)^{3/2}e^{B_{\mathrm{D}}/T}.\tag{14}
$$

Both the neutron and proton density are proportional to the baryon density, so we have \(n_n \simeq n_p \simeq n_b = \eta n_\gamma^{(0)}\). Substituting this into Eq. (14), we obtain:

$$
\frac{n_{\mathrm{D}}}{n_{\mathrm{b}}}\sim\eta_{\mathrm{b}}\left(\frac{T}{m_{p}}\right)^{3/2}e^{B_{\mathrm{D}}/T}.\tag{15}
$$

As long as \(B_D/T\) is not too large, the small prefactor dominates this expression.

### 4.2.1 Neutron abundance

We begin by solving for the neutron–proton ratio. Protons can be converted into neutrons via weak interactions, \(p + e^- \rightarrow n + \nu_e\) for example.
From Eq. (7), the proton/neutron equilibrium ratio in the nonrelativistic limit (so that \(E_i (p) = m_i + p^2/2m_i\) ) is

$$
\frac{n_p^{(0)}}{n_n^{(0)}}=\frac{e^{-m_p/T}\int dp\mathrm{~}p^2e^{-p^2/2m_pT}}{e^{-m_n/T}\int dp\mathrm{~}p^2e^{-p^2/2m_nT}}.\tag{16}
$$

We know that the integrals in Eq.(16) are proportional to \(m^{3/2}\), so we can neglect the mass difference as the ratio \(m_p/m_n\) approaches 1. However, we can neglect the difference in the exponents, and we are left with:

$$
\frac{n_p^{(0)}}{n_n^{(0)}}=e^{\mathcal{Q}/T}\tag{17}
$$

with \(\mathcal{Q} \equiv m_n - m_p = 1.293MeV\). Therefore, at high temperatures, there are as many neutrons as protons. As the temperature drops beneath 1 MeV, the neutron fraction goes down. If weak interactions operated efficiently enough to maintain equilibrium indefinitely, then it would drop to zero.
Define:

$$
X_n \equiv \frac{n_n}{n_n + n_p},\tag{18}
$$

that is, \(X_n\) is the ratio of neutrons to total nuclei.In equilibrium,

$$
X_n\to X_{n,\mathrm{EQ}}\equiv\frac{1}{1+n_p^{(0)}/n_n^{(0)}}.\tag{19}
$$

To track the evolution of \(X_n\), let us start from Eq. (10), with 1 = neutron, 3 = proton, and 2,  4 = leptons in complete equilibrium (\(n_l = n_l^{(0)}\)). Then,

$$
a^{-3}\frac{d(n_na^3)}{dt}=n_l^{(0)}\langle\sigma v\rangle\left\{\frac{n_pn_n^{(0)}}{n_p^{(0)}}-n_n\right\}.\tag{20}
$$

We have already determined the ratio $n_n^{( 0) }/ n_p^{( 0) }= e^{- Q/ T}$ and we can identify $n_l^{( 0) }\langle \sigma v\rangle$ as $\lambda _{np}$, the rate for neutron$\to$proton conversion since it multiplies $n_n$ in the loss term. Also, if we rewrite $n_n$ on the left as $(n_n+n_p)X_n$, then the total density times $a^3$ can be taken outside the derivative(for baryon number conservation), leaving

$$
\frac{dX_n}{dt}=\lambda_{np}\left\{(1-X_n)e^{-\mathcal{Q}/T}-X_n\right\}.\tag{21}
$$

Define:

$$
x \equiv \frac{\mathcal{Q}}{T},\tag{22}
$$

the left of Eq.(21) becomes $\dot{x} d X_n/dx$, so we need an expression for $\dot{x} = -x \dot{T}/T$. Since $T \propto a^{-1}$, we have

$$
\frac{1}{T}\frac{dT}{dt}=-H=-\sqrt{\frac{8\pi G\rho}{3}},\tag{23}
$$

Nucleosynthesis occurs in the radiationdominated era, so the main contribution to the energy density $\rho$ comes from relativistic particles. Recall from Ch. 2 that the contribution to the energy density from relativistic particles is

$$
\begin{aligned}
\rho & =\frac{\pi^2}{30}T^4\left[\sum_{s=\mathrm{bosons}}g_s+\frac{7}{8}\sum_{s=\mathrm{fermions}}g_s\right]\quad\text{(s relativistic)} \\
 & \equiv g_*\frac{\pi^2}{30}T^4.
\end{aligned}\tag{24}
$$

The effective numbers of relativistic degrees of freedom, $g_*$, is a function of the temperature. At temperatures of order l MeV, the contributing species are: photons $(g_\gamma=2)$, neutrinos $(g_\nu=6)$,and electrons and positrons $(g_{e^+}=g_{e^-}=2).$ Adding up leads to $g_*\simeq10.75$, roughly constant throughout the regime of interest. Then, Eq. (21) becomes

$$
\frac{dX_n}{dx}=\frac{x\lambda_{np}}{H(x=1)}\left\{e^{-x}-X_n(1+e^{-x})\right\}.\tag{25}
$$

### 4.2.2 Light elements abundances

To simplify the understanding and calculation, we think light element production occurs instantaneously at a temperature $T_{nuc}$ when the energetics compensates for the small baryon-to-photon ratio, which means that the temperature is so cold that the energy of almost photons is smaller than the binding energy of this light element. Let us consider deuterium as an example, with Eq.(15):

$$ln(\eta_b) + \frac{3}{2} ln(T_{nuc}/m_p) \sim -\frac{B_D}{T_{nuc}}.\tag{26}$$

Helium is produced almost immediately after deuterium, so almost all remaining neutrons at $T \sim T_{nuc}$ are processed into $^4He$. Since each helium nucleus has 2 protons and 2 neutrons, the final helium abundance is equal to half the neutron abundance at $T_{nuc}$. Therefore, we have:

$$Y_p \equiv \frac{4n_{He}}{n_b} \simeq 2X_n(T_{nuc}),\tag{27}$$
which yields a final helium mass fraction of 0.22. One important feature of this result is that it depends only weakly on the baryon-to-photon ratio $\eta_b$. This is because the neutron fraction $X_n$ depends only logarithmically on $\eta_b$ through $T_{nuc}$. Not all of the deuterium produced at $T_{nuc}$ is converted into helium, but only a small fraction of it. A trace amount remains unburned, simply because the reaction that eliminates it, $D + p \rightarrow ^3He + \gamma$,is not completely efficient. While deuterium is depleted via these reactions after $T_\mathrm{nuc}$, it eventually freezes out at a mass fraction of order $\bar{3}\times10^{-5}.$ If the baryon density is low, then the reactions proceed more slowly, and the depletion is not as effective. Therefore, low baryon density inevitably results in more deuterium; the sensitivity is quite stark. This fact, combined with the possibility of measuring deuterium in high-redshift gas clouds by looking for absorption in the spectra of distant QSOs (see Sect. l.3), turns the deuterium abundance into an important probe of the baryon density.
## 4.3 Recombination
