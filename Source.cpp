#include <iostream>
#include <conio.h>
#include <iomanip>
#include <cmath>
#include <limits>
#include <string>
#include <CImg.h>

using namespace cimg_library;
using namespace std;



#define pi 3.141592

#define byte unsigned char
#define ush  unsigned short
#define uint unsigned int
#define ll   long long
#define ull  unsigned long long

#define Image CImg<unsigned char>

#define XSize 300
#define YSize 300
#define Magnification 2



Image       main_image(XSize * Magnification, YSize * Magnification, 1, 3);
CImgDisplay main_display(main_image, "Simulation");



template <typename datatype>
struct LPart {
	LPart<datatype>* pointer;
	datatype payload;
};
template <typename datatype>
class list {
private:
	LPart<datatype>* start;
	bool empty = 1;

public:
	int length = 0;

	list();
	void append(datatype dlt);
	void add(datatype dlt); // it`s like append, but adds new element to the start of the list (instead of end)
	void del(short ind);
	datatype pop();
	datatype* getlink(short ind);
	std::string str(std::string sep);

	datatype& operator[](short ind) {
		LPart<datatype>* last = start;
		for (short q = 0; q < ind; q++) {
			last = (*last).pointer;
		}

		return (*last).payload;
	}
};
template <typename datatype>
list<datatype>::list() {

}
template <typename datatype>
void list<datatype>::append(datatype dlt) {
	if (empty) {
		start = new LPart<datatype>;
		(*start).payload = dlt;
		(*start).pointer = start;

		empty = 0;
	}
	else {
		if (length == 0) {
			(*start).payload = dlt;
			(*start).pointer = start;
		}
		else {
			LPart<datatype>* appendix = new LPart<datatype>;
			(*appendix).payload = dlt;
			(*appendix).pointer = appendix;

			LPart<datatype>* last = start;
			while ((*last).pointer != last) {
				last = (*last).pointer;
			}
			(*last).pointer = appendix;
		}
	}

	length += 1;
}
template <typename datatype>
void list<datatype>::add(datatype dlt) {
	if (empty) {
		start = new LPart<datatype>;
		(*start).payload = dlt;
		(*start).pointer = start;

		empty = 0;
	}
	else {
		if (length == 0) {
			(*start).payload = dlt;
			(*start).pointer = start;
		}
		else {
			LPart<datatype>* appendix = new LPart<datatype>;
			(*appendix).payload = dlt;
			(*appendix).pointer = start;

			start = appendix;
		}
	}

	length += 1;
}
template <typename datatype>
void list<datatype>::del(short ind) {
	length -= 1;

	if (ind == 0) {
		start = start->pointer;
	}
	else {
		LPart<datatype>* part = start;
		LPart<datatype> newpointer;

		for (short q = 0; q < ind - 1; q++) {
			part = (*part).pointer;
		}

		(*part).pointer = (*(*part).pointer).pointer;
	}
}
template <typename datatype>
datatype list<datatype>::pop() {
	length -= 1;

	datatype out = (*start).payload;
	start = (*start).pointer;
	return out;
}
template <typename datatype>
datatype* list<datatype>::getlink(short ind) {
	LPart<datatype>* last = start;
	for (short q = 0; q < ind; q++) {
		last = (*last).pointer;
	}

	return &(*last).payload;
}
template <typename datatype>
std::string list<datatype>::str(std::string sep) {
	string out = "";

	if (length > 0) {
		for (short q = 0; q < length; q++) {
			out += get(q);
			out += sep;
		}
		out.erase(out.end() - 1);
	}

	return out;
}

struct vector {
	double x = 0;
	double y = 0;
};
class phasor {
public:
	double magnitude;
	double phase;

	phasor() {
		magnitude = 0;
		phase = 0;
	}
	phasor(vector v) {
		magnitude = sqrt(v.x * v.x + v.y * v.y);
		phase = atan(v.x / v.y);
	}
	void operator=(phasor other) {
		this->magnitude = other.magnitude;
		this->phase = other.phase;
	}
};



double output_power_coefficient = 1;
double output_dissipation_coefficient = 1;
double output_power_to_dissipation_ratio_coefficient = 1;
double output_power_on_target_distance_coefficient = 1;
double output_dissipation_on_target_distance_coefficient = 1;
double output_power_to_dissipation_ratio_on_target_distance_coefficient = 1;

bool ShowInLogarithmicScale = 1;

class Shape {
private:
	double getintensity(double x1, double y1, double z1) {
		double x_val = 0;
		double y_val = 0;

		double x2, y2;  // z2 is always zero
		double range_squear;
		double phase_shift;

		double r;
		double A = 0;
		double dA = 0.1;
		while (A < 2 * pi) {
			double dr = 1;

			r = 0;
			while (r < d1 / 2) {
				x2 = r * cos(A) + XSize / 2;
				y2 = r * sin(A) + YSize / 2;

				range_squear = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + z1 * z1;
				phase_shift = 2 * pi * sqrt(range_squear) / L;

				x_val += r * dA * dr * cos(phase_shift) / range_squear;
				y_val += r * dA * dr * sin(phase_shift) / range_squear;

				r += dr;
			}

			r = R - d2 / 2;
			while (r < R + d2 / 2) {
				x2 = r * cos(A) + XSize / 2;
				y2 = r * sin(A) + YSize / 2;

				range_squear = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + z1 * z1;
				phase_shift = 2 * pi * sqrt(range_squear) / L;

				x_val += r * dA * dr * cos(phase_shift) / range_squear;
				y_val += r * dA * dr * sin(phase_shift) / range_squear;

				r += dr;
			}

			A += dA;
		}

		return sqrt(x_val * x_val + y_val * y_val);
	}
	phasor getphasor(double x1, double y1, double z1) {
		vector out;
		out.x = 0;
		out.y = 0;

		double x2, y2;  // z2 is always zero
		double range_squear;
		double phase_shift;

		double r;
		double A = 0;
		double dA = 0.1;
		while (A < 2 * pi) {
			double dr = 1;

			r = 0;
			while (r < d1 / 2) {
				x2 = r * cos(A) + XSize / 2;
				y2 = r * sin(A) + YSize / 2;

				range_squear = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + z1 * z1;
				phase_shift = 2 * pi * sqrt(range_squear) / L;

				out.x += r * dA * dr * cos(phase_shift) / range_squear;
				out.y += r * dA * dr * sin(phase_shift) / range_squear;

				r += dr;
			}

			r = R - d2 / 2;
			while (r < R + d2 / 2) {
				x2 = r * cos(A) + XSize / 2;
				y2 = r * sin(A) + YSize / 2;

				range_squear = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + z1 * z1;
				phase_shift = 2 * pi * sqrt(range_squear) / L;

				out.x += r * dA * dr * cos(phase_shift) / range_squear;
				out.y += r * dA * dr * sin(phase_shift) / range_squear;

				r += dr;
			}

			A += dA;
		}

		return phasor(out);
	}
	double getintensityafterwall(double x1, double y1, double d1, phasor* intensity_list) {
		double x_val = 0;
		double y_val = 0;

		phasor p;
		double range_squear;
		double xcor, ycor;
		double dA = 0.1;
		for (double A = 0; A < 2 * pi; A += dA) {
			for (ush r = 0; r < hole_radius; r++) {
				p = intensity_list[r];

				xcor = r * cos(A) + XSize / 2;
				ycor = r * sin(A) + YSize / 2;

				if (p.magnitude != 0) {
					range_squear = (xcor - x1) * (xcor - x1) + (ycor - y1) * (ycor - y1) + d1 * d1;

					p.phase += 2 * pi * sqrt(range_squear) / L;

					x_val += dA * r * sin(p.phase) * p.magnitude / range_squear;
					y_val += dA * r * cos(p.phase) * p.magnitude / range_squear;
				}
			}
		}

		return sqrt(x_val * x_val + y_val * y_val);
	}
	phasor getphasorafterwall(double x1, double y1, double d1, phasor* intensity_list) {
		vector out;
		out.x = 0;
		out.y = 0;

		phasor p;
		double range_squear;
		double xcor, ycor;
		double dA = 0.1;
		for (double A = 0; A < 2 * pi; A += dA) {
			for (ush r = 0; r < hole_radius; r++) {
				p = intensity_list[r];

				xcor = r * cos(A) + XSize / 2;
				ycor = r * sin(A) + YSize / 2;

				if (p.magnitude != 0) {
					range_squear = (xcor - x1) * (xcor - x1) + (ycor - y1) * (ycor - y1) + d1 * d1;

					p.phase += 2 * pi * sqrt(range_squear) / L;

					out.x += dA * r * sin(p.phase) * p.magnitude / range_squear;
					out.y += dA * r * cos(p.phase) * p.magnitude / range_squear;
				}
			}
		}

		return phasor(out);
	}
	vector getintensityonwall(double x1, double y1) {
		vector out;
		out.x = 0;
		out.y = 0;

		double x2, y2;  // z2 is always zero
		double range_squear;
		double phase_shift;

		double r;
		double A = 0;
		double dA = 0.1;
		while (A < 2 * pi) {
			double dr = 1;

			r = 0;
			while (r < d1 / 2) {
				x2 = r * cos(A) + XSize / 2;
				y2 = r * sin(A) + YSize / 2;

				range_squear = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + H * H;
				phase_shift = 2 * pi * sqrt(range_squear) / L;

				out.x += r * dA * dr * cos(phase_shift) / range_squear;
				out.y += r * dA * dr * sin(phase_shift) / range_squear;

				r += dr;
			}

			r = R - d2 / 2;
			while (r < R + d2 / 2) {
				x2 = r * cos(A) + XSize / 2;
				y2 = r * sin(A) + YSize / 2;

				range_squear = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + H * H;
				phase_shift = 2 * pi * sqrt(range_squear) / L;

				out.x += r * dA * dr * cos(phase_shift) / range_squear;
				out.y += r * dA * dr * sin(phase_shift) / range_squear;

				r += dr;
			}

			A += dA;
		}

		return out;
	}
	double calc_goodness() {
		if (check_validity()) {
			for (ush r = 0; r < hole_radius; r++) {
				intensity_profile[r] = getintensityonwall(XSize / 2, YSize / 2 + r);
			}
			output_power = 0;
			output_dissipation = 0;
			for (double d = 0; d < L * 20; d += L / 4) {
				output_dissipation += pow(getintensityafterwall(XSize / 2 + L * 4, YSize / 2, d, intensity_profile) / 0.02, 3);
			}
			for (double d = 2 * L; d < L * 20; d += L / 4) {
				output_power += pow(getintensityafterwall(XSize / 2, YSize / 2, d, intensity_profile) / 0.02, 3);
			}
			output_power_on_target_distance = 0;
			output_dissipation_on_target_distance = 0;
			double dA = 0.05;
			double power_area_angle = pi / 18;
			double target_distance = L * 20;
			for (double A = 0; A < power_area_angle; A += dA) {
				output_power_on_target_distance += target_distance * target_distance * A * dA * pow(getintensityafterwall(XSize / 2,
					YSize / 2 + target_distance * sin(A),
					target_distance * cos(A),
					intensity_profile) / 0.02, 3);
			}
			for (double A = power_area_angle; A < pi / 2; A += dA) {
				output_dissipation_on_target_distance += target_distance * target_distance * A * dA * pow(getintensityafterwall(XSize / 2,
					YSize / 2 + target_distance * sin(A),
					target_distance * cos(A),
					intensity_profile) / 0.02, 3);
			}

			goodness = 0;
			goodness += output_power * 0.03 * output_power_coefficient;
			goodness += output_power_on_target_distance * 20000 * output_power_on_target_distance_coefficient;
			goodness -= output_dissipation * 1840 * output_dissipation_coefficient;
			goodness -= output_dissipation_on_target_distance *  32800 * output_dissipation_on_target_distance_coefficient;
			goodness += 0.025 * output_power / output_dissipation * output_power_to_dissipation_ratio_coefficient;
			goodness += 2840 * pow(output_power_on_target_distance / output_dissipation_on_target_distance, 2) * output_power_to_dissipation_ratio_on_target_distance_coefficient;
		}
		else {
			goodness = -INT32_MAX;
		}

		return goodness;
	}
	double* getval(byte ind) {
		if (ind == 0) {
			return &R;
		}
		if (ind == 1) {
			return &H;
		}
		if (ind == 2) {
			return &d1;
		}
		if (ind == 3) {
			return &d2;
		}
		if (ind == 4) {
			return &hole_radius;
		}
	}

	phasor intensity_profile[XSize * 2];

public:
	double R = 26.46;              // ring hole radius
	double d1 = 37.5;              // diameter of the inner hole
	double d2 = 3.616;             // width of the ring hole
	double L = 15;                 // wavelength of the soundwave
	double H = 6.548;              // viewing distance
	double hole_radius = 45;       // output hole radius
						   /*
						   there is a flat sound absorbing surface
						   inside that surface there is a circular hole with diameter d1
						   also there is a ring - shaped hole with inner radius R - d2 / 2 and outer radius R + d2 / 2
						   the inner hole is positioned in the center of the ring hole
						   also there is a point sound source
						   */

	double goodness;
	double output_power;
	double output_dissipation;
	double output_power_on_target_distance;
	double output_dissipation_on_target_distance;

	Shape() {

	}
	void print_data() {
		calc_goodness();

		cout << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl;
		cout << "Viewing height:                                            " << H << endl;
		cout << endl;
		cout << "Inner hole diameter:                                       " << d1 << endl;
		cout << "Ring hole radius:                                          " << R << endl;
		cout << "Ring hole width:                                           " << d2 << endl;
		cout << "Output hole diameter:                                      " << hole_radius * 2 << endl;
		cout << "Wavelength:                                                " << L << endl;
		cout << endl;
		cout << "Relative goodness:                                         " << goodness << endl;
		cout << "Wave after wall power:                                     " << output_power << endl;
		cout << "Wave after wall dissipaton:                                " << output_dissipation << endl;
		cout << "Wave after wall power on target distance:                  " << output_power_on_target_distance << endl;
		cout << "Wave after wall dissipaton on target distance:             " << output_dissipation_on_target_distance << endl;
		cout << "Wave power to dissipation ratio:                           " << output_power / output_dissipation << endl;
		cout << "Wave power to dissipation ratio on target distance:        " << output_power_on_target_distance / output_dissipation_on_target_distance << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl;
		cout << endl;
	}
	void operator=(Shape other) {
		this->R = other.R;
		this->H = other.H;
		this->L = other.L;
		this->d1 = other.d1;
		this->d2 = other.d2;
		this->hole_radius = other.hole_radius;
	}
	void find_local_maximum(ush number_of_iterations, double optimisation_value, bool is_independent, double opimisation_value_delta) {
		if (is_independent) {
			cout << "*********************************************************" << endl;
			cout << "Start data:" << endl;
			print_data();
			cout << endl << endl;
		}

		bool change;
		for (ush iteration = 0; iteration < number_of_iterations; iteration++) {
			change = 0;

			Shape s;
			for (byte value_id = 0; value_id < 5; value_id++) {
				for (char c = 0; c < 1; c++) {
					s = *this;
					s.getval(value_id);
					(*s.getval(value_id)) *= 1 + optimisation_value * (c * 2 - 1);

					s.calc_goodness();
					if (s.goodness > this->goodness) {
						*this = s;
						change = 1;

						calc_goodness();
						cout << "New goodness: " << goodness << endl;
					}
				}
			}

			if (!change) {
				if (is_independent) {
					optimisation_value /= 1 + optimisation_value * opimisation_value_delta;
				}
				else {
					break;
				}
			}
			if (optimisation_value == 0) {
				if (is_independent) {
					cout << endl << "+++++++++++++++++++++++++++++ Maximum found! +++++++++++++++++++++++++++++" << endl << endl;
				}
				break;
			}
		}

		if (is_independent) {
			cout << endl << "End optimisation value: " << optimisation_value << endl;

			cout << endl << endl;
			cout << "End data:" << endl;
			print_data();
			cout << "*********************************************************" << endl;
		}
	}
	void successive_approach(double start_optimisation_value, double end_optimisation_value, double step) {
		cout << endl << "Optimisation started" << endl << endl;

		cout << "Start data:" << endl;
		print_data();

		for (double optimisation_coefficient = start_optimisation_value; optimisation_coefficient > end_optimisation_value; optimisation_coefficient -= step) {
			find_local_maximum(25, optimisation_coefficient, 0, 0);
		}

		cout << "End data:" << endl;
		print_data();

		cout << endl << "Optimisation ended" << endl << endl;
	}
	void get_intensity_in_point() {
		cout << "Wait, until the image shows up" << endl;
		show_after_wall();

		cout << "Tap on the picture to get intensity in that point" << endl;
		cout << "Press 'S' to stop" << endl;

		bool brk = 0;
		while (1) {
			if (GetKeyState('S') & 0x8000) {
				break;
			}
			else {
				while (!main_display.button() & 1) {
					Sleep(5);
					if (GetKeyState('S') & 0x8000) {
						brk = 1;
						break;
					}
				}
				if (brk) {
					break;
				}

				double x = main_display.mouse_x() / Magnification;
				double y = main_display.mouse_y() / Magnification;
				cout << "Intensity in chosen point " << "(" << x << " " << y << ")" << ": " << getintensityafterwall(XSize / 2, y, x, intensity_profile) << endl;
				Sleep(100);
			}
		}
	}
	void picking_method() {
		Shape other;

		double Rstart = L; double Rend = L * 3; double dR = L / 3;
		double Hstart = L / 2; double Hend = L * 2; double dH = L / 4;
		double d1start = L / 2; double d1end = L * 3; double dd1 = L / 2;
		double d2start = L / 15; double d2end = L / 3; double dd2 = L / 15;
		double nhrstart = L / 2; double nhrend = L * 4; double dnhr = L / 2;

		calc_goodness();
		for (double new_R = Rstart; new_R < Rend; new_R += L / 2) {
			for (double new_H = Hstart; new_H < Hend; new_H += L) {
				for (double new_d1 = d1start; new_d1 < d1end; new_d1 += L / 2) {
					for (double new_d2 = d2start; new_d2 < d2end; new_d2 += L / 10) {
						for (double new_hole_radius = nhrstart; new_hole_radius < nhrend; new_hole_radius += dnhr) {
							other.R = new_R;
							other.H = new_H;
							other.d1 = new_d1;
							other.d2 = new_d2;
							other.hole_radius = new_hole_radius;

							cout << "Attempt data: " << new_R << ' ' << new_H << ' ' << new_d1 << ' ' << new_d2 << ' ' << new_hole_radius << endl;

							other.calc_goodness();
							if (other.goodness > goodness) {
								cout << "Found new best!" << endl;

								*this = other;
								successive_approach(0.1, 0.01, 0.01);
								show_after_wall();
							}
						}
					}
				}
			}
		}
	}
	bool check_validity() {
		if (H < L / 2 - 1) {
			return 0;
		}
		if (R - d2 / 2 <= d1 / 2) {
			return 0;
		}
		if (hole_radius <= L / 2) {
			return 0;
		}
		if (hole_radius >= XSize * 2) {
			return 0;
		}

		return 1;
	}
	double show_slice() {
		cout << "Loading image..." << endl;

		byte color[3];
		double intensity;
		for (ush xcor = 0; xcor < XSize; xcor++) {
			for (ush ycor = 0; ycor < YSize; ycor++) {
				intensity = getintensity(xcor, ycor, H);

				if ((((abs(xcor - XSize / 2) == int(hole_radius)) || (abs(ycor - YSize / 2) == int(hole_radius))) && ((ycor == YSize / 2) || xcor == XSize / 2))) {
					color[0] = 0;
					color[1] = 0;
					color[2] = 0;
				}
				else {
					double val;
					if (ShowInLogarithmicScale) {
						val = log10(intensity * 3000 + 10);
					}
					else {
						val = pow(intensity, 0.15) * pi;
					}
					color[0] = 127 * (1 + sin(val * 0.5 + pi * 0));
					color[1] = 127 * (1 + sin(val * 1 + pi * 0));
					color[2] = 117 * (1 + sin(val * 0.5 + pi * 0.5));
				}

				for (ush dx = 0; dx < Magnification; dx++) {
					for (ush dy = 0; dy < Magnification; dy++) {
						main_image.draw_point(xcor * Magnification + dx, ycor * Magnification + dy, color, 1);
					}
				}
			}
		}
		main_image.display(main_display);

		calc_goodness();

		return goodness;
	}
	double show_slice_phase() {
		cout << "Loading image..." << endl;

		byte color[3];
		phasor p;
		for (ush xcor = 0; xcor < XSize; xcor++) {
			for (ush ycor = 0; ycor < YSize; ycor++) {
				if ((((abs(xcor - XSize / 2) == int(hole_radius)) || (abs(ycor - YSize / 2) == int(hole_radius))) && ((ycor == YSize / 2) || xcor == XSize / 2))) {
					color[0] = 0;
					color[1] = 0;
					color[2] = 0;
				}
				else {
					p = getphasor(xcor, ycor, H);

					color[0] = 127 * (1 + sin(p.phase));
					color[1] = 127 * (1 + sin(p.phase + pi / 2));
					color[2] = 117 * (1 + sin(p.phase + pi));
				}

				for (ush dx = 0; dx < Magnification; dx++) {
					for (ush dy = 0; dy < Magnification; dy++) {
						main_image.draw_point(xcor * Magnification + dx, ycor * Magnification + dy, color, 1);
					}
				}
			}
		}
		main_image.display(main_display);

		calc_goodness();

		return goodness;
	}
	double show_profile() {
		cout << "Loading image..." << endl;

		byte color[3];
		double intensity;
		for (ush h = 1; h < XSize + 1; h++) {
			for (ush ycor = 0; ycor <= YSize / 2; ycor++) {
				intensity = getintensity(XSize / 2, ycor, h);

				if ((h == int(H)) & ((ycor == YSize / 2 - int(hole_radius)) | (ycor == YSize / 2 + int(hole_radius)))) {
					color[0] = 0;
					color[1] = 0;
					color[2] = 0;
				}
				else {
					double val;
					if (ShowInLogarithmicScale) {
						val = log10(intensity * 3000 + 10);
					}
					else {
						val = pow(intensity, 0.15) * pi;
					}
					color[0] = 127 * (1 + sin(val * 0.5 + pi * 0));
					color[1] = 127 * (1 + sin(val * 1 + pi * 0));
					color[2] = 117 * (1 + sin(val * 0.5 + pi * 0.5));
				}

				for (ush dx = 0; dx < Magnification; dx++) {
					for (ush dy = 0; dy < Magnification; dy++) {
						main_image.draw_point((h - 1) * Magnification + dx, ycor * Magnification + dy, color, 1);
						main_image.draw_point((h - 1) * Magnification + dx, YSize * Magnification - ycor * Magnification - dy, color, 1);
					}
				}
			}
		}
		main_image.display(main_display);

		calc_goodness();

		return goodness;
	}
	double show_profile_phase() {
		cout << "Loading image..." << endl;

		byte color[3];
		phasor p;
		for (ush h = 1; h < XSize + 1; h++) {
			for (ush ycor = 0; ycor <= YSize / 2; ycor++) {
				if ((h == int(H)) & ((ycor == YSize / 2 - int(hole_radius)) | (ycor == YSize / 2 + int(hole_radius)))) {
					color[0] = 0;
					color[1] = 0;
					color[2] = 0;
				}
				else {
					p = getphasor(XSize / 2, ycor, h);

					color[0] = 127 * (1 + sin(p.phase));
					color[1] = 127 * (1 + sin(p.phase + pi / 2));
					color[2] = 117 * (1 + sin(p.phase + pi));
				}

				for (ush dx = 0; dx < Magnification; dx++) {
					for (ush dy = 0; dy < Magnification; dy++) {
						main_image.draw_point((h - 1) * Magnification + dx, ycor * Magnification + dy, color, 1);
						main_image.draw_point((h - 1) * Magnification + dx, YSize * Magnification - ycor * Magnification - dy, color, 1);
					}
				}
			}
		}
		main_image.display(main_display);

		calc_goodness();

		return goodness;
	}
	double show_after_wall() {
		cout << "Loading image..." << endl;

		for (ush r = 0; r < hole_radius; r++) {
			intensity_profile[r] = getintensityonwall(XSize / 2, YSize / 2 + r);
		}

		byte color[3];
		double intensity;
		for (ush d = 1; d < XSize + 1; d++) {
			for (ush ycor = 0; ycor <= YSize / 2; ycor++) {
				if ((d == 2) & ((ycor == YSize / 2 - int(hole_radius)) | (ycor == YSize / 2 + int(hole_radius)))) {
					color[0] = 0;
					color[1] = 0;
					color[2] = 0;
				}
				else {
					intensity = getintensityafterwall(XSize / 2, ycor, d, intensity_profile);

					double val;
					if (ShowInLogarithmicScale) {
						val = log10(intensity * 3000 + 10);
					}
					else {
						val = pow(intensity, 0.15) * pi;
					}
					color[0] = 127 * (1 + sin(val * 0.5 + pi * 0));
					color[1] = 127 * (1 + sin(val * 1 + pi * 0));
					color[2] = 117 * (1 + sin(val * 0.5 + pi * 0.5));
				}

				for (ush dx = 0; dx < Magnification; dx++) {
					for (ush dy = 0; dy < Magnification; dy++) {
						main_image.draw_point((d - 1) * Magnification + dx, ycor * Magnification + dy, color, 1);
						main_image.draw_point((d - 1) * Magnification + dx, YSize * Magnification - ycor * Magnification - dy, color, 1);
					}
				}
			}
		}
		main_image.display(main_display);

		calc_goodness();
		print_data();

		return 0;
	}
	double show_after_wall_phase() {
		cout << "Loading image..." << endl;

		for (ush r = 0; r < hole_radius; r++) {
			intensity_profile[r] = getintensityonwall(XSize / 2, YSize / 2 + r);
		}

		byte color[3];
		phasor p;
		for (ush d = 1; d < XSize + 1; d++) {
			for (ush ycor = 0; ycor <= YSize / 2; ycor++) {
				if ((d == 2) & ((ycor == YSize / 2 - int(hole_radius)) | (ycor == YSize / 2 + int(hole_radius)))) {
					color[0] = 0;
					color[1] = 0;
					color[2] = 0;
				}
				else {
					p = getphasorafterwall(XSize / 2, ycor, d, intensity_profile);

					color[0] = 127 * (1 + sin(p.phase));
					color[1] = 127 * (1 + sin(p.phase + pi / 2));
					color[2] = 117 * (1 + sin(p.phase + pi));
				}

				for (ush dx = 0; dx < Magnification; dx++) {
					for (ush dy = 0; dy < Magnification; dy++) {
						main_image.draw_point((d - 1) * Magnification + dx, ycor * Magnification + dy, color, 1);
						main_image.draw_point((d - 1) * Magnification + dx, YSize * Magnification - ycor * Magnification - dy, color, 1);
					}
				}
			}
		}
		main_image.display(main_display);

		calc_goodness();
		print_data();

		return 0;
	}
};
Shape main_shape;



void Showsettings() {
	cout << "Current settings:" << endl;
	cout << "    Image settings:" << endl;
	cout << "	      0.0 show in logarithmic scale: " << ShowInLogarithmicScale << endl;
	cout << "              // if 1, wave simulations will be shown in logarithmic scale. Else, they will be shown in pow(0.15) scale." << endl;
	cout << endl;
	cout << "    Goodness calculation settings:" << endl;
	cout << "        1.0 Output average power coefficient:                                                    " << output_power_coefficient << endl;
	cout << "        1.1 Output average dissipation coefficient:                                              " << output_dissipation_coefficient << endl;
	cout << "        1.2 Output average power to average dissipation ratio coefficient:                       " << output_power_to_dissipation_ratio_coefficient << endl;
	cout << "        1.3 Output power on target distance coefficient:                                         " << output_power_coefficient << endl;
	cout << "        1.4 Output dissipation on target distance coefficient:                                   " << output_dissipation_coefficient << endl;
	cout << "        1.5 Output power on target distance to dissipation on target distance ratio coefficient: " << output_power_to_dissipation_ratio_coefficient << endl;
	cout << "             // goodness is a value, that determines, how good the shape is." << endl;
	cout << "             // it`s a key component of artificial evolution that you use in this program" << endl;
	cout << "             // to calculate goodness, the program calculates shape properties, for example average output power" << endl;
	cout << "             // goodness is a sum of various shape properties, multiplied by an appropriate coefficient" << endl;
	cout << "             // so, if you change coefficients, what you do is you change evolution criteria." << endl;

}
void Settings() {
	Showsettings();

	cout << endl;
	cout << "To exit settings, enter 'exit'" << endl;
	cout << "To change a program parameter, enter parameter code." << endl;
	cout << "For example, 1.4 means 'Output dissipation on target distance  coefficient'" << endl;
	cout << endl;

	string command;
	while (1) {
		cout << "Enter number: ";
		getline(cin, command);

		if (command == "exit") {
			break;
		}
		else if (command == "0.0") {
			cin >> ShowInLogarithmicScale;
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin.clear();
		}
		else if (command == "1.0") {
			cin >> output_power_coefficient;
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin.clear();
		}
		else if (command == "1.1") {
			cin >> output_dissipation_coefficient;
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin.clear();
		}
		else if (command == "1.2") {
			cin >> output_power_to_dissipation_ratio_coefficient;
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin.clear();
		}
		else if (command == "1.3") {
			cin >> output_power_on_target_distance_coefficient;
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin.clear();
		}
		else if (command == "1.4") {
			cin >> output_dissipation_on_target_distance_coefficient;
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin.clear();
		}
		else if (command == "1.5") {
			cin >> output_power_to_dissipation_ratio_on_target_distance_coefficient;
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin.clear();
		}
		else {
			cout << "Wrong number" << endl;
		}
	}
}
void Showhelp() {
	cout << endl;
	cout << "This program is created to calculate a desirable properties of a device," << endl;
	cout << "that i call ring concelntrator. It consists of a sound source of wavelength L" << endl;
	cout << "inside a sound absorbing box." << endl;
	cout << "There are two holes in that box. One of them is a ring - shaped hole. Radius of this hole is R" << endl;
	cout << "and the width of this hole is d2. In the center of the ring hole there is another hole." << endl;
	cout << "It`s circle - shaped and has a diameter d1. Together this two holes creare a wave, that" << endl;
	cout << "interferes with itself, creating areas, where the wave completely cancels itself out." << endl;
	cout << "If you now place a barier with a circular hole with right diameter, theoretically, you" << endl;
	cout << "shoud get a direct soundwave. The distance from the hole to the barier is H." << endl;
	cout << endl;
	cout << "List of commands:" << endl;
	cout << "    help" << endl;
	cout << "        shows help" << endl << endl;
	cout << "    exit" << endl;
	cout << "        breaks the program" << endl << endl;
	cout << "    settings" << endl;
	cout << "        show and change program settings" << endl << endl;
	cout << "    save image" << endl;
	cout << "        saves current screen image in a file" << endl << endl;
	cout << "    get intensity in point" << endl;
	cout << "        you tap on a point on a simulation and get intensity in that point" << endl << endl;
	cout << "    show data" << endl;
	cout << "        shows shape data" << endl << endl;
	cout << "    show inner wave slice" << endl;
	cout << "        shows wave slice on chosen H" << endl << endl;
	cout << "    show inner wave slice phase" << endl;
	cout << "        shows wave phase slice on chosen H" << endl << endl;
	cout << "    show inner wave profile" << endl;
	cout << "        shows wave profile" << endl << endl;
	cout << "    show inner wave profile phase" << endl;
	cout << "        shows wave phase profile before the wall" << endl << endl;
	cout << "    show output wave profile" << endl;
	cout << "        shows wave profile after exiting the wall" << endl << endl;
	cout << "    show output wave profile phase" << endl;
	cout << "        shows wave phase profile after exiting the wall" << endl << endl;
	cout << "    successive approach" << endl;
	cout << "        relatively fast and effective way to optimize the shape to a closest maximum" << endl << endl;
	cout << "    picking method" << endl;
	cout << "        slowly, but surely finds best shape properties" << endl << endl;
	cout << "    local maximum" << endl;
	cout << "        a very fast way to bring your shape to the closest maximum. not reccomended to use, since it does not work very well." << endl << endl;
	cout << "    H" << endl;
	cout << "        changes value of viewing height" << endl << endl;
	cout << "    R" << endl;
	cout << "        changes value of ring hole radius" << endl << endl;
	cout << "    r" << endl;
	cout << "        changes value of output hole radius" << endl << endl;
	cout << "    d1" << endl;
	cout << "        changes value of circle hole diameter" << endl << endl;
	cout << "    d2" << endl;
	cout << "        changes value of ring hole width" << endl << endl;
	cout << "    L" << endl;
	cout << "        changes value of wavelength" << endl;
	cout << endl;
	cout << "Red color on wave visualisation means, that the amplitude of the wave is high" << endl;
	cout << "Green color on wave visualisation means, that the amplitude of the wave is medium" << endl;
	cout << "Red color on wave visualisation means, that the amplitude of the wave is low" << endl;
	cout << "Black dots on visualisations represent wall edges" << endl;
	cout << endl;
}
void Execute(std::string command) {
	string input_buffer;

	if (command == "help") {
		Showhelp();
	}
	else if (command == "exit") {
		exit(0);
	}
	else if (command == "settings") {
		Settings();
	}
	else if (command == "save image") {
		main_image.save("Simulation.bmp", time(0));
	}
	else if (command == "get intensity in point") {
		main_shape.get_intensity_in_point();
	}
	else if (command == "show data") {
		main_shape.print_data();
	}
	else if (command == "show inner wave slice") {
		main_shape.show_slice();
		main_shape.print_data();
	}
	else if (command == "show inner wave slice phase") {
		main_shape.show_slice_phase();
		main_shape.print_data();
	}
	else if (command == "show inner wave profile") {
		main_shape.show_profile();
		main_shape.print_data();
	}
	else if (command == "show inner wave profile phase") {
		main_shape.show_profile_phase();
		main_shape.print_data();
	}
	else if (command == "show output wave profile") {
		main_shape.show_after_wall();
	}
	else if (command == "show output wave profile phase") {
		main_shape.show_after_wall_phase();
	}
	else if (command == "local maximum") {
		ush number_of_iterations;
		double optimisation_value;
		cout << "Enter number of iterations: ";
		cin >> number_of_iterations;
		cin.ignore(numeric_limits<streamsize>::max(), '\n');
		cin.clear();
		cout << "Enter the optimisation value: ";
		cin >> optimisation_value;
		cin.ignore(numeric_limits<streamsize>::max(), '\n');
		cin.clear();
		main_shape.find_local_maximum(number_of_iterations, optimisation_value, 1, 0.3);
	}
	else if (command == "successive approach") {
		cout << "'Local maximum' function tries to optimize the shape using a given parameter, named 'optimisation value'." << endl;
		cout << "It determines, how much will the shape change during evolution process. High OV will optimize the shape" << endl;
		cout << "very fast, but roughly. Low OV will work slowly, but accurately. To ensure both quality and speed, this" << endl;
		cout << "function uses various optimisation value. First it`s hight (START), then, during the optimisation" << endl;
		cout << "process, every iteration it gets a bit lower (STEP), until it reaches the minimal value (END). Then the" << endl;
		cout << "cycle stops." << endl << endl;

		double start, end, step;
		while (1) {
			cout << "Enter START:    " << endl;
			cin >> start;
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin.clear();
			cout << "Enter END:      " << endl;
			cin >> end;
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin.clear();
			cout << "Enter STEP:     " << endl;
			cin >> step;
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin.clear();

			if (step > 0 & start > end) {
				break;
			}
			else {
				cout << "You entered wrong values. STEP must be higher than zero. START must be higher than END." << endl;
			}
		}
		main_shape.successive_approach(start, end, step);
	}
	else if (command == "picking method") {
		main_shape.picking_method();
	}
	else if (command == "H") {
		cout << endl;
		cout << "Enter new H: ";
		cin >> main_shape.H;
		cin.ignore(numeric_limits<streamsize>::max(), '\n');
		cin.clear();
	}
	else if (command == "R") {
		cout << endl;
		cout << "Enter new R: ";
		cin >> main_shape.R;
		cin.ignore(numeric_limits<streamsize>::max(), '\n');
		cin.clear();
	}
	else if (command == "r") {
		cout << endl;
		cout << "Enter new r: ";
		cin >> main_shape.hole_radius;
		cin.ignore(numeric_limits<streamsize>::max(), '\n');
		cin.clear();
	}
	else if (command == "d1") {
		cout << endl;
		cout << "Enter new d1: ";
		cin >> main_shape.d1;
		cin.ignore(numeric_limits<streamsize>::max(), '\n');
		cin.clear();
	}
	else if (command == "d2") {
		cout << endl;
		cout << "Enter new d2: ";
		cin >> main_shape.d2;
		cin.ignore(numeric_limits<streamsize>::max(), '\n');
		cin.clear();
	}
	else if (command == "L") {
		cout << endl;
		cout << "Enter new L: ";
		cin >> main_shape.L;
		cin.ignore(numeric_limits<streamsize>::max(), '\n');
		cin.clear();
	}
	else {
		cout << endl << "Wrong command. Enter 'help' for help." << endl;
		return;
	}

	cout << endl << "Execution went successful" << endl;
}



int main() {
	main_image.fill(0);

	cout << std::setprecision(24);

	cout << "Enter 'help' for help" << endl;

	string command;
	while (1) {
		cout << ">>> ";
		getline(cin, command);

		Execute(command);
	}
}