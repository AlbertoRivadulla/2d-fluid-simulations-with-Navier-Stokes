# 2-dimensional fluid simulations with the Navier-Stokes equation

Simulations of 2-dimensional fluids solving the Navier-Stokes equations.

### To do
- `BoundaryConditions`
    - List of points with boundary conditions, for the obstacles.
        - Each of these will be a struct with the indices of the point, and the
        direction: N, S, E, W, NE, SE, NW and SW
        - The method `applyObstaclesBCs` will iterate over these and apply them
        one by one
    - Implement the method to add a rectangle obstacle.
    - Finish the method to add a circular obstacle.
- `NavierStokesSolver`
    - In the method `solvePoisson`, when computing the initial value of the
    pressure, use the correct value for the number of fluid cells.

