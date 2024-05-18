// Your First C++ Program

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

class Body {
public:
    float x_vel_new = 0;
    float y_vel_new = 0;
    float x_pos_new = 0;
    float y_pos_new = 0;
    float x_pos;
    float y_pos;
    float x_vel;
    float y_vel;
    float mass;
    float x_acc = 0;
    float y_acc = 0;
    // Initialise Body Object
    Body(float x_pos_in, float y_pos_in, float x_vel_in, float y_vel_in, float mass_in):
        x_pos(x_pos_in), y_pos(y_pos_in), x_vel(x_vel_in), y_vel(y_vel_in), mass(mass_in) {}
};
void write_out(std::ofstream& out, std::vector<Body> &bodies, float time) {
    for (int i = 0; i < bodies.size(); ++i){
        out << i << ','
        << time << ','
        << bodies[i].x_pos << ','
        << bodies[i].y_pos << ','
        << bodies[i].x_vel << ','
        << bodies[i].y_vel << ',' << '\n';
    }
}

float calc_dst_magnitude(float xi, float yi, float xj, float yj) {
    return pow(pow((xi - xj), 2) + pow((yi - yj), 2), 0.5);
}

float calc_acc(float G, float dst, float pos1, float pos2, float m2) {
    return G * m2 * (pos2 - pos1) / pow(dst, 3);
}

std::vector<float> evalk(float G, std::vector<Body> &bodies,
    int curr_idx, float b1_xpos, float b1_ypos, float b1_xvel, float b1_yvel, float delta_t) {
    float acc_x = 0;
    float acc_y = 0;
    for (int j = 0; j < bodies.size(); ++j) {
        if (curr_idx == j) {
            continue;
        }
        Body b2 = bodies[j];
        float dst = calc_dst_magnitude(b1_xpos, b1_ypos, b2.x_pos, b2.y_pos);
        acc_x += calc_acc(G, dst, b1_xpos, b2.x_pos, b2.mass);
        acc_y += calc_acc(G, dst, b1_ypos, b2.y_pos, b2.mass);
    }
    return {b1_xvel, b1_yvel, acc_x, acc_y};
}

void rk4(std::vector<Body> &bodies, float delta_t, float G) {
  // Calc loop
  std::vector<std::vector<std::vector<float>>> ks;
  for (int i = 0; i < bodies.size(); ++i) {
    Body b1 = bodies[i];
    //b1.x_pos_new = b1.x_pos;
   // b1.y_pos_new = b1.y_pos;
    //b1.x_vel_new = b1.x_vel;
    //b1.y_vel_new = b1.y_vel;
    std::vector<float> k1 = evalk(G, bodies, i, b1.x_pos, b1.y_pos, b1.x_vel, b1.y_vel, delta_t);
    std::vector<float> k2 = evalk(G, bodies, i, b1.x_pos + 0.5 * k1[0], b1.y_pos + 0.5 * k1[1], b1.x_vel + 0.5 * k1[2], b1.y_vel + 0.5 * k1[3], delta_t);
    std::vector<float> k3 = evalk(G, bodies, i, b1.x_pos - 0.5*k1[0] + 0.5 * k2[0], b1.y_pos -0.5*k1[1] + 0.5 * k2[1], b1.x_vel -0.5*k1[2]+ 0.5 * k2[2], b1.y_vel -0.5*k1[3]+ 0.5 * k2[3], delta_t);
    std::vector<float> k4 = evalk(G, bodies, i, b1.x_pos - 0.5*k2[0]+ k3[0], b1.y_pos -0.5*k2[1]+ k3[1], b1.x_vel -0.5*k2[2] + k3[2], b1.y_vel -0.5*k2[2]+ k3[3], delta_t);
    std::vector<float> krenew = evalk(G, bodies, i, b1.x_pos -k3[0], b1.y_pos -k3[1], b1.x_vel - k3[2], b1.y_vel -k3[3], delta_t);
    ks.push_back({k1, k2, k3, k4});
    // Updating loop

  }
  for (int i = 0; i < ks.size(); ++i){
    for (int j = 0; j < bodies.size(); ++j) {
        bodies[j].x_pos += delta_t * ((ks[i][0][0] + ks[i][3][0]) / 6.0 + (ks[i][1][0] + ks[i][2][0]) / 3.0);
        bodies[j].y_pos += delta_t * ((ks[i][0][1] + ks[i][3][1]) / 6.0 + (ks[i][1][1] + ks[i][2][1]) / 3.0);
        bodies[j].x_vel += delta_t * ((ks[i][0][2] + ks[i][3][2]) / 6.0 + (ks[i][1][2] + ks[i][2][2]) / 3.0);
        bodies[j].y_vel += delta_t * ((ks[i][0][3] + ks[i][3][3]) / 6.0 + (ks[i][1][3] + ks[i][2][3]) / 3.0);
    }}
  //for (int i = 0; i < bodies.size(); ++i) {
  //      bodies[i].x_pos = bodies[i].x_pos_new;
  //      bodies[i].y_pos = bodies[i].y_pos_new;
  //      bodies[i].x_vel = bodies[i].x_vel_new;
  //       bodies[i].y_vel = bodies[i].y_vel_new;
  //  }
}

int main() {
    /** Assume parameters.txt exists */
    std::ifstream input_file("parameters.txt");
    std::ofstream output_file("output.csv");
    std::string firstline;
    std::string G;
    std::string T;
    std::string delta_t;
    getline(input_file, G, ' ');
    getline(input_file, T, ' ');
    getline(input_file, delta_t);
    std::vector<Body> bodies;

    std::string xpos_i;
    std::string ypos_i;
    std::string xvel_i;
    std::string yvel_i;
    std::string m_i;
    std::vector<float> x_vel_track;
    std::vector<float> y_vel_track;
    while (!input_file.eof()) {
        getline(input_file,xpos_i, ' ');
        getline(input_file,ypos_i, ' ');
        getline(input_file,xvel_i, ' ');
        getline(input_file,yvel_i, ' ');
        getline(input_file,m_i);
        Body new_body = Body(stof(xpos_i), stof(ypos_i), stof(xvel_i), stof(yvel_i), stof(m_i));
        bodies.push_back(new_body);
    }
    // Do your RK4 here
    for (float t = 0; t < stof(T); t += stof(delta_t)) {
        rk4(bodies, stof(delta_t), stof(G));
        write_out(output_file, bodies, t);
    }
    // for (unsigned int i = 0; i < xpos_i_l.size(); i++) {
    //     write_out(output_file,
    //         xpos_i_l[i],
    //         ypos_i_l[i],
    //         xvel_i_l[i],
    //         yvel_i_l[i],
    //         m_i_l[i]);
    // }
    return 0;
}
