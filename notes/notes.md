## SSVF algorithm pseudo-code
#### Step 1: vertex finding
```python
# List of coupled tracks (two-tracks vertex approximation)
couples = []
indices = [i for i in range(len(tracks))]
# While there are tracks to be coupled
while len(tracks)>1:
    # Select the track to be coupled
    i,t = indices.pop(0),tracks.pop(0)
    #selecting the closest track
    j = argmin(distances)
    # if tracks are close enough, couple them
    if distances[j] < threshold:
        tc = tracks.pop(j)
        indices.pop(j)
        couples.append([t,tc])
return couples
```

#### Step 2: vertex fitting
```python
# fitting loop
while True:
    # Find the SV that minimizes the chi2
    # from the current selected tracks
    SVfitted = min(chi2(tracks))
    # If chi2 is good enough, break
    if chi2 < threshold:
        break
    # Else, reject the worst track
    tracks.pop(argmax(chi2s))
    # If there are less than 2 tracks,
    # no vertex has been found
    if len(tracks)<2:
        return None
return SVfitted
```

### Tracks distances calculations

Let the track line be represented in a parametric form:
$$
r: \Bigg \{ \begin{array}{l}
    x = x_0 + \cos(\varphi)\sin(\theta)  t \\
    y = y_0 + \sin(\varphi)\sin(\theta) t \\
    z = z_0 + \cos(\theta)t \\
    \end{array}
$$
where $\bm{P}=(x_0,y_0,z_0)^{\top}$ is the origin of the track and $$\bm{\hat{v}} = (\cos(\varphi)\sin(\theta),\sin(\varphi)\sin(\theta),\cos(\theta))^{\top}$$
is the track versor.
#### Distance between two tracks
The minimum distance between two lines is the length of the only segment that is orthogonal to both of them. This segment is the solution of the following system:
$$\Bigg \{ \begin{array}{l}
    \Delta\bm{P} \cdot \hat{r}_1 = 0 \\
    \Delta\bm{P} \cdot \hat{r}_2 = 0 \\
\end{array}$$
where $\Delta \bm{P}$ is the vectorial difference between the parametric points of the two lines, which are:
$$      
Pr_1 = (x_0 + at, y_0 + bt, z_0 +ct)\\
Pr_2 = (x^{'}_0 + a^{'}s, y^{'}_0 + b^{'}s, z^{'}_0 +c^{'}s)\\
$$
Thus the system can be written as:
$$
\Bigg \{ \begin{array}{l}
    \Delta \bm{P}_0 \cdot \bm{\hat{r}}_1 + (\bm{\hat{r}}_1 \cdot \bm{\hat{r}}_2)s - \bm{\hat{r}}^{2}_1t = 0\\
    \Delta \bm{P}_0 \cdot \bm{\hat{r}}_2 + \bm{\hat{r}}^{2}_2s - (\bm{\hat{r}}_1 \cdot \bm{\hat{r}}_2)t = 0\\            
\end{array}
$$
Where $\Delta \bm{P}_0$ is the vector that connects the tracks' origins.
The system is solved by:
$$
\Bigg \{ \begin{array}{l}
    t^* = -\frac{1}{\sin^2(\alpha)}\big[\cos(\alpha)(\Delta \bm{P}_0 \cdot \bm{\hat{r}}_2) - \Delta \bm{P}_0 \cdot \bm{\hat{r}}_1 \big] \\
    s^* = \cos(\alpha)t^* - \Delta \bm{P}_0 \cdot \bm{\hat{r}}_2
\end{array}
$$
By plugging the values of $t^*$ and $s^*$ in the parametric representation of the lines one finds the lines' closest points. The distance between the two points is the minimum distance between the lines.

#### Distance between a track and a point
The distance between a line and a point $\bm{P}_d=(x_d,y_d,z_d)^{\top}$ can be evaluated in two steps. The first step is  finding the plane orthogonal to the line that contains the point.
A general plane is given by:
$$\pi: ax+by+cz+d=0$$
A plane orthogonal to the line must have its versor parallel to the line's one: $\bm{\hat{\pi}} = c \bm{\hat{v}}$. WLOG, $c=1$.
So a general plane orthogonal to the line is given by:
$$v_0 x + v_1 y + v_2 z + d = 0$$
Now, the plane must contain $\bm{P}_d$ so:
$$d^* = -(v_0 x_d + v_1 y_d + v_2 z_d)$$

The line's point that's closest to $\bm{P}_d$ is the intersection between the plane and the line. By substituting the parametric representation of the line in the plane equation:
$$
v_0 x_0 + v_1 y_0 + v_2 z_0 + v_0^2 t + v_1^2 t + v_2^2 t + d^* = 0
$$
Finally:
$$t^* = -(\frac{v_0x_0 + v_1 y_0 + v_0z_2}{v_0^2 + v_1^2 + v_2^2} + d^*)$$

The distance is the distance between $\bm{P}_d$ and the parametric representation of the line evaluated at $t = t^*$.