# SPH Fluid Simulation — Navier-Stokes in C++

2D Smoothed Particle Hydrodynamics (SPH) simulation in real time, built in C++ with [raylib](https://www.raylib.com/). The particle physics include actual Navier-Stokes terms with their proper kernels (density, pressure, viscosity).

The Navier-Stokes equation solution is one of the [Millennium Prize Problems](https://en.wikipedia.org/wiki/Millennium_Prize_Problems), to which anyone who solves it will receive one million US dollars. The Navier-Stokes equation has no current general solution but is currently the most accurate way of illustrating the motion and physics behind fluids.

There is a toggle between honey and water to illustrate differences in the Navier-Stokes terms (specifically viscosity, as honey is extremely viscous). This toggle happens in real time by pressing either "w" for water or "h" for honey.

There is also attraction and repulsion at the mouse pointer for interaction with the fluid (you can generate your own waves or droplets).

## Features

- **1000 particles** simulated at 60-80 fps
- **Full SPH physics**: density estimation, pressure via equation of state, pressure forces, and viscosity — each derived from a relative kernel
- **Spatial hashing** — a uniform grid reduces the neighbour search from O(n²) to roughly O(n), recovering most of the FPS lost to the physics calculations. This produced nearly a double boost in FPS.
- **Live fluid switching** between water (low viscosity, thin) and honey (high viscosity, thick) in real time
- **Speed-based colour mapping** — blue-white for water, orange-yellow for honey. This gives an idea of the particle movement.
- **Mouse interaction** — left-click repels particles, right-click attracts them

## Build

### Prerequisites
- [raylib](https://github.com/raysan5/raylib) installed on your system

### Windows (MinGW)
```bash
g++ fluid_sim.cpp -o fluid_sim.exe -lraylib -lopengl32 -lgdi32 -lwinmm -std=c++17
fluid_sim.exe
```

### Linux / macOS
```bash
g++ fluid_sim.cpp -o fluid_sim -lraylib -lm -lpthread -ldl -std=c++17
./fluid_sim
```

## Controls

| Input | Effect |
|---|---|
| `w` | Toggle Water |
| `h` | Toggle Honey |
| Left Mouse Button | Repel force from mouse pointer |
| Right Mouse Button | Attract force from mouse pointer |
| Close Window / Esc | Quit |

## The physics

The simulation involves a simplified form of the Navier-Stokes equation (one of the Millennium Prize Problems yet to be solved — our current best illustration of the motion of particles in fluids):

```
∂u/∂t = -(1/ρ)∇p + ν∇²u + g
```
ρ is density, p is pressure, ν is dynamic viscosity (resistance to flow), g is the external force (gravity, in this case), and u is the velocity vector of the fluid.

Each frame, every particle computes:

**1. Density** — a Poly6-style smoothing kernel sums mass contributions from all neighbours within the interaction radius:
```
ρᵢ = Σⱼ mⱼ · (1 - r/h)²
```

**2. Pressure** — an equation of state converts local density excess into pressure, clamped to never go negative (fluids resist compression, not expansion):
```
pᵢ = k · max(ρᵢ - ρ₀, 0)
```

**3. Pressure force** — a Spiky-style gradient kernel pushes particles apart proportionally to combined pressure, preventing particles from over-compressing.

**4. Viscosity** — a Laplacian-style kernel averages velocity with nearby particles, producing the syrupy motion that distinguishes water from honey. Nearby particles "share" each other's velocity evenly.

## Performance: spatial hashing

A naive implementation checks every particle against every other particle each frame — O(n²), or 3 million pair checks per frame at 1000 particles across three physics passes. This implementation instead buckets particles into a uniform grid sized to the interaction radius, and each particle only checks the 3×3 block of cells around it, nearly doubling FPS.

### A debugging note worth sharing

After adding spatial hashing, there was barely any boost in FPS, even though checks per frame should have dropped dramatically. Going to 3000 particles made FPS collapse rather than stay flat — a clear sign the grid wasn't actually limiting the neighbour search. Digging into this surfaced three separate issues:

1. **Non-uniform particle distribution.** Resulting from the particle physics, fluid pools settle at the bottom due to pressure from particles in the middle of the fluid. While the average particles per cell stayed around 4-5 (what we want), the busiest cells at the bottom had around 20-30 particles per cell, which meant most of the actual neighbour checks were still happening in those overloaded cells — undermining the benefit of the hashing exactly where it mattered most.

2. **A density/rest_density scale mismatch.** The Poly6-style kernel produces small density values (around 4-10), but `rest_density` was initially set far above that range. This meant the equation of state never triggered — pressure stayed near 0 even when particles were crushed right next to each other. Recalibrating `rest_density` to match the kernel's actual output scale fixed this: pressure activated correctly, particles spread to a more realistic resting distance, and the spatial hashing grid finally delivered its expected FPS boost.

3. **An exponentially small viscosity kernel.** After getting the simulation running well by adjusting other constants, testing extreme viscosity values produced no visible change at all. The cause was a `powf(interaction_radius, 6)` term in the original Laplacian-style kernel — a normalisation constant designed for physically accurate SPH at molecular scale, where you're dealing with enormous numbers of tiny particles. At pixel scale with only ~1000 particles, that term made the viscosity force evaluate to roughly 10⁻¹⁰ — effectively zero, no matter how high the viscosity coefficient was set. The fix was to replace it with a simplified linear falloff, `(1 - r/h)`, which produces sensible values at this scale. Viscosity worked correctly after that.

This is a good example of how a simulation can look visually correct while quietly hiding incorrect physics underneath.

## Potential extensions

- 3D SPH with a proper renderer (instanced spheres or marching cubes for a surface mesh)
- Surface tension term
- Variable particle mass/size for multi-fluid mixing
- GPU compute shader for the neighbour search at much higher particle counts
