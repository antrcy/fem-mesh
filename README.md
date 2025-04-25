

On souhaite implémenter la fonction qui fait la transformation affine entre deux triangles : au lieu de travailler en base canonique, on travaille avec les polynômes de Lagrange. On suppose $\hat x_i$ les sommets de l'élément de référence, envoyés sur l'élément réel $K$ de sommets $x_i$. La transformation géométrique est de la forme :

$$
T_K(\hat x) = \sum _{i = 1}^3 x_i \hat \phi_i (\hat x)
$$

On a donc :

$$
T_K(\hat x_j) = x_j
$$

L'application affine $\hat x \to B_K \hat x + c_k$ telle que :

$$
T_K(\hat K) = K
$$

On calcule sa jacobienne pour obtenir :

$$
D T_K = 
\begin{pmatrix}
    -x_1^1 + x_2^1 & -x_1^1 + x_3^1 \\
    -x_1^2 + x_2^2 & -x_1^2 + x_3^2
\end{pmatrix}
$$

Avec $x_i = (x_i^1, x_i^2)$. Vu qu'elle est affine, on a :

$$
B_K = DT_K
$$

et on identifie :

$$
b = x_1
$$

On a donc une formule explicite, qui n'est que combinaison linéaire des fonctions de base sur l'élément de référence :

$$
\hat \phi_1 = 1 - x - y
$$

$$
\hat \phi_2 = x
$$

$$
\hat \phi_3 = y
$$

Les gradients se calculent directement :

$$
\nabla \hat \phi_1 = \begin{pmatrix} -1 \\ -1 \end{pmatrix}
$$

$$
\nabla \hat \phi_2 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}
$$

$$
\nabla \hat \phi_3 = \begin{pmatrix} 0 \\ 1 \end{pmatrix}
$$

## Calcul des normes

L'objectif est de transformer toutes nos intégrations sur un élément réel en une intégration sur un élément de référence.

### Calcul de la norme L2

On calcule l'erreur L2 entre la solution et son approximation sur l'élément $K$ :

$$
\int_K \big( \text{sol}(x, y) - f(x, y) \big)^2 dxdy
$$

$$
 = \int_K \big( \sum_{i = 1}^3 u_i^K \phi_i^K(x, y) - f(x, y) \big)^2 dxdy
$$

$$
 = \int_{\hat K} \big( \sum_{i = 1}^3 u_i^K \hat \phi _i (x, y) - f(T_K (x, y ))\big)^2 |\det (\nabla T_K)| dxdy
$$

### Calcul de la norme H1

Le terme manquant est calculer comme suit :

$$
\int_K \| \nabla \text{sol}(x, y) - \nabla f(x, y) \|^2 dxdy
$$

$$
= \int_K \| \sum_{i = 1}^3 u_i^K \nabla \phi_i^K(x, y) - \nabla f(x, y) \|^2 dxdy
$$

$$
= \int_{\hat K} \| \sum_{i = 1}^3 u_i ^K \nabla \phi_i^K(T_K(x, y)) - \nabla f(T_K(x, y)) \|^2 |\det(\nabla T_K)| dxdy
$$

Le terme $\nabla f(T_K(x, y))$ n'est pas développé. En revanche, $\nabla \phi_i^K(T_K(x, y))$ n'est pas entièrement connu, puisqu'on ne connaît que $\hat \phi_i$. On utilise la proposition suivante :

$$
\nabla \phi_i^K(T_K(x, y)) = \nabla \hat \phi_i \nabla T_K(x, y)^{-T} 
$$

Tous les termes sont connus. On obtient finalement :

$$
\int_{\hat K} \| \sum_{i = 1}^3 u_i ^K \nabla \hat \phi_i \nabla T_K(x, y)^{-T}  - \nabla f(T_K(x, y)) \|^2 |\det(\nabla T_K)| dxdy
$$
