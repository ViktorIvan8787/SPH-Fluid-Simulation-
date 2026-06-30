#include "raylib.h"
#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <string>

// Create for each particle
struct Particle {
    Vector2 position;
    Vector2 velocity;
    float radius;
    float density;
    float pressure;
};


int main() {

    // Simulation window attributes
    const int windowWidth = 1000;
    const int windowHeight = 750;

    InitWindow(windowWidth, windowHeight, "Fluid Simulation");
    SetTargetFPS(80);

    // =================== Simulation constants ====================
    int number_of_particles = 1000;

    float gravity = 500.0f;
    float bounce_coefficient = 1.0f;
    float drag_coefficient = 0.99f;
    float interaction_radius = 60.0f;

    float mouse_radius = 80.0f;
    float mouse_strength = 10000.0f;

    float mass = 1.0f;
    float viscosity, stiffness, rest_density;

    // Input choice window between honey and water
    // fluid_type = 0 for water, fluid_type = 1 for honey
    bool fluid_chosen = false;
    int fluid_type = 0;
    std::string fluid;

    while (!fluid_chosen && !WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(BLACK);

        DrawText("Choose a fluid:", 300, 300, 30, WHITE);
        DrawText("Press W to simulate Water", 300, 350, 24, WHITE);
        DrawText("Press H to simulate Honey", 300, 390, 24, WHITE);
        DrawText("Navier-Stokes Implementation", 790, 10, 14, WHITE);
        DrawText("(Density, Pressure, and Velocity from SPH kernels)", 640, 28, 14, WHITE);
        DrawText("Interaction with Spatial Hashing", 780, 48, 14, WHITE);

        if (IsKeyPressed(KEY_W)) {
            fluid_type = 0;
            fluid_chosen = true;
        }
        if (IsKeyPressed(KEY_H)) {
            fluid_type = 1;
            fluid_chosen = true;
        }

        EndDrawing();
    }

    if (fluid_type == 0) {
        // Water
        viscosity = 0.05f;
        stiffness = 1000.0f;
        rest_density = 10.0f;
    } else {
        // Honey
        viscosity = 2.5f;
        stiffness = 1000.0f;
        rest_density = 10.0f;
    }

    // Spatial hashing constants
    int size_cell = (int)interaction_radius;
    int grid_columns = windowWidth / size_cell + 1;
    int grid_rows = windowHeight / size_cell + 1;

    std::vector<std::vector<int>> grid(grid_columns * grid_rows);

    // Creating many particles
    std::vector<Particle> particles;

    for (int i = 0; i < number_of_particles; i++) {
        Particle p;
        p.position = { (float)GetRandomValue(50, 950), (float)GetRandomValue(400, 780) };
        p.velocity = { (float)GetRandomValue(-50, 50), 0.0f };
        p.radius = 4.0f;
        particles.push_back(p);
    }

    // ==================== Simulation Loop ===================
    while (!WindowShouldClose()) {

        // Allow fluid type to be changed live
        if (IsKeyPressed(KEY_W)) fluid_type = 0;
        if (IsKeyPressed(KEY_H)) fluid_type = 1;

        if (fluid_type == 0) {
            viscosity = 0.05f;
            stiffness = 1000.0f;
            rest_density = 10.0f;
        } else {
            viscosity = 2.5f;
            stiffness = 1000.0f;
            rest_density = 10.0f;
        }

        float dt = GetFrameTime();

        // Gravity + integration + boundaries
        for (int i = 0; i < number_of_particles; i++) {
            Particle& p = particles[i];

            p.velocity.y += gravity * dt;
            p.velocity.y *= drag_coefficient;
            p.velocity.x *= drag_coefficient;

            p.position.x += p.velocity.x * dt;
            p.position.y += p.velocity.y * dt;

            if (p.position.y + p.radius > windowHeight) {
                p.position.y = windowHeight - p.radius - 2;
                p.velocity.y *= -bounce_coefficient;
            }
            if (p.position.y - p.radius < 0.0f) {
                p.position.y = p.radius + 2;
                p.velocity.y *= -bounce_coefficient;
            }
            if (p.position.x - p.radius < 0.0f) {
                p.position.x = p.radius + 2;
                p.velocity.x *= -bounce_coefficient;
            }
            if (p.position.x + p.radius > windowWidth) {
                p.position.x = windowWidth - p.radius - 2;
                p.velocity.x *= -bounce_coefficient;
            }
        }

        // ===== Spatial hashing grid =====
        for (auto& cell : grid) cell.clear();

        for (int i = 0; i < number_of_particles; i++) {
            int x_cell = (int)(particles[i].position.x / size_cell);
            int y_cell = (int)(particles[i].position.y / size_cell);
            int cell_number = y_cell * grid_columns + x_cell;
            grid[cell_number].push_back(i);
        }

        // === Density (Poly6 kernel) ===
        for (int i = 0; i < number_of_particles; i++) {
            particles[i].density = 0.0f;

            int x_cell = (int)(particles[i].position.x / size_cell);
            int y_cell = (int)(particles[i].position.y / size_cell);

            for (int temp_x_cell = x_cell - 1; temp_x_cell <= x_cell + 1; temp_x_cell++) {
                for (int temp_y_cell = y_cell - 1; temp_y_cell <= y_cell + 1; temp_y_cell++) {

                    if (temp_x_cell < 0 || temp_x_cell >= grid_columns || temp_y_cell < 0 || temp_y_cell >= grid_rows) continue;

                    int cell_number = temp_y_cell * grid_columns + temp_x_cell;

                    for (int j : grid[cell_number]) {
                        float distance_x = particles[j].position.x - particles[i].position.x;
                        float distance_y = particles[j].position.y - particles[i].position.y;
                        float dist = sqrtf(distance_x * distance_x + distance_y * distance_y);

                        if (dist < interaction_radius) {
                            float W = powf(1.0f - dist / interaction_radius, 2);
                            particles[i].density += mass * W;
                        }
                    }
                }
            }
        }

        // === Pressure (equation of state) ===
        for (int i = 0; i < number_of_particles; i++) {
            particles[i].pressure = stiffness * std::max(particles[i].density - rest_density, 0.0f);
        }

        // === Pressure forces (simplified Spiky gradient) ===
        for (int i = 0; i < number_of_particles; i++) {
            int x_cell = (int)(particles[i].position.x / size_cell);
            int y_cell = (int)(particles[i].position.y / size_cell);

            for (int temp_x = x_cell - 1; temp_x <= x_cell + 1; temp_x++) {
                for (int temp_y = y_cell - 1; temp_y <= y_cell + 1; temp_y++) {

                    if (temp_x < 0 || temp_x >= grid_columns || temp_y < 0 || temp_y >= grid_rows) continue;

                    int cell_number = temp_y * grid_columns + temp_x;

                    for (int j : grid[cell_number]) {
                        if (i == j) continue;

                        Particle& a = particles[i];
                        Particle& b = particles[j];

                        Vector2 delta = { b.position.x - a.position.x, b.position.y - a.position.y };
                        float distance_two_p = sqrt(delta.x * delta.x + delta.y * delta.y);

                        if (distance_two_p < interaction_radius && distance_two_p > 0) {
                            float gradient_pressure = -2.0f * (1.0f - distance_two_p / interaction_radius) / interaction_radius;

                            float vx = delta.x / distance_two_p;
                            float vy = delta.y / distance_two_p;

                            float force = mass * (particles[i].pressure + particles[j].pressure) * gradient_pressure;

                            particles[i].velocity.x += force * vx * dt;
                            particles[i].velocity.y += force * vy * dt;
                        }
                    }
                }
            }
        }

        // === Viscosity (simplified Laplacian kernel) ===
        for (int i = 0; i < number_of_particles; i++) {
            int x_cell = (int)(particles[i].position.x / size_cell);
            int y_cell = (int)(particles[i].position.y / size_cell);

            for (int temp_x = x_cell - 1; temp_x <= x_cell + 1; temp_x++) {
                for (int temp_y = y_cell - 1; temp_y <= y_cell + 1; temp_y++) {

                    if (temp_x < 0 || temp_x >= grid_columns || temp_y < 0 || temp_y >= grid_rows) continue;

                    int cell_number = temp_y * grid_columns + temp_x;

                    for (int j : grid[cell_number]) {
                        if (i == j) continue;

                        Particle& a = particles[i];
                        Particle& b = particles[j];

                        Vector2 delta = { b.position.x - a.position.x, b.position.y - a.position.y };
                        float distance_two_p = sqrt(delta.x * delta.x + delta.y * delta.y);

                        if (distance_two_p < interaction_radius && distance_two_p > 0) {
                            float lapW = (1.0f - distance_two_p / interaction_radius);

                            particles[i].velocity.x += viscosity * mass * (b.velocity.x - a.velocity.x) * lapW * dt;
                            particles[i].velocity.y += viscosity * mass * (b.velocity.y - a.velocity.y) * lapW * dt;
                        }
                    }
                }
            }
        }

        // === Mouse interaction ===
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            Vector2 mouse = GetMousePosition();
            for (int i = 0; i < number_of_particles; i++) {
                Particle& p = particles[i];
                Vector2 delta = { p.position.x - mouse.x, p.position.y - mouse.y };
                float distance_from_mouse = sqrt(delta.x * delta.x + delta.y * delta.y);

                if (distance_from_mouse < mouse_radius && distance_from_mouse > 0.0f) {
                    Vector2 direction = { delta.x / distance_from_mouse, delta.y / distance_from_mouse };
                    float closeness = 1.0f - (distance_from_mouse / mouse_radius);
                    float force = closeness * mouse_strength;

                    p.velocity.x += direction.x * force * dt;
                    p.velocity.y += direction.y * force * dt;
                }
            }
        }

        if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
            Vector2 mouse = GetMousePosition();
            for (int i = 0; i < number_of_particles; i++) {
                Particle& p = particles[i];
                Vector2 delta = { p.position.x - mouse.x, p.position.y - mouse.y };
                float distance_from_mouse = sqrt(delta.x * delta.x + delta.y * delta.y);

                if (distance_from_mouse < mouse_radius && distance_from_mouse > 0.0f) {
                    Vector2 direction = { delta.x / distance_from_mouse, delta.y / distance_from_mouse };
                    float closeness = 1.0f - (distance_from_mouse / mouse_radius);
                    float force = closeness * mouse_strength;

                    p.velocity.x += -direction.x * force * dt;
                    p.velocity.y += -direction.y * force * dt;
                }
            }
        }

        // === Drawing ===
        BeginDrawing();
        ClearBackground(BLACK);

        for (int i = 0; i < number_of_particles; i++) {
            float speed = sqrtf(particles[i].velocity.x * particles[i].velocity.x + particles[i].velocity.y * particles[i].velocity.y);
            float t = std::min(speed / 300.0f, 1.0f);

            Color col;
            if (fluid_type == 0) {
                // Water: deep blue -> white
                col = {
                    (unsigned char)(0 + t * 255),
                    (unsigned char)(20 + t * 235),
                    (unsigned char)(150 + t * 105),
                    220
                };
            } else {
                // Honey: dark orange -> bright yellow
                col = {
                    (unsigned char)(180 + t * 75),
                    (unsigned char)(90 + t * 145),
                    (unsigned char)(0 + t * 20),
                    220
                };
            }

            DrawCircleV(particles[i].position, particles[i].radius, col);
        }

        DrawFPS(10, 10);

        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT) || IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
            DrawCircleLinesV(GetMousePosition(), mouse_radius, WHITE);
        }

        fluid = (fluid_type == 0) ? "Water" : "Honey";

        DrawText(TextFormat("Current fluid type: %s (press h or w to change)", fluid.c_str()), 10, 35, 18, WHITE);
        DrawText("Right click to attract, Left click to repel", 10, 55, 18, WHITE);
        DrawText("Navier-Stokes Implementation", 790, 10, 14, WHITE);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}
