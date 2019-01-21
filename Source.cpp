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
CImgDisplay main_display(main_image);



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



double contrast_coefficient = 1;
double size_allignment_coefficient = 1;
double size_constancy_coefficient = 1;
double output_power_coefficient = 1;
double output_dissipation_coefficient = 1;

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

		return sqrt(x_val * x_val + y_val * y_val);  // return value from 0 to 1
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

		return sqrt(x_val * x_val + y_val * y_val);  // return value from 0 to 1
	}
	vector getintensityonwall(double x1, double y1) {
		vector out;

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
		get_radius(H, 1);

		if (check_validity()) {
			wave_r1 = get_radius(H * 2, 0);
			contrast = 1 - pow(minimal_intensity / maximal_intensity, 2);
			size_allignment = 1 - abs(hole_radius - R) / (hole_radius + R);
			size_constancy = pow(min(hole_radius, wave_r1) / max(hole_radius, wave_r1), 2);

			for (ush r = 0; r < hole_radius; r++) {
				intensity_profile[r] = getintensityonwall(XSize / 2, YSize / 2 + r);
			}
			output_power = 0;
			output_dissipation = 0;
			for (ush d = 0; d < L * 100; d += L / 4) {
				output_power += pow(getintensityafterwall(XSize / 2, YSize / 2, L * 5, intensity_profile) / 0.02, 3);
				output_dissipation += pow(getintensityafterwall(XSize / 2 + hole_radius * 2, YSize / 2, L * 5, intensity_profile) / 0.02, 3);
			}

			goodness = 0;
			goodness += contrast * 3 * contrast_coefficient;
			goodness += size_allignment * 1 * size_allignment_coefficient;
			goodness += size_constancy * 10 * size_constancy_coefficient;
			goodness += output_power * 0.01 * output_power_coefficient;
			goodness -= output_dissipation * 500 * output_dissipation_coefficient;
		}
		else {
			hole_radius = 0;
			wave_r1 = 0;
			goodness = 0;
			maximal_intensity = 0;
			minimal_intensity = 0;
			contrast = 0;
			size_allignment = 0;
			size_constancy = 0;
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
	double R = 13.7777;            // ring hole radius
	double d1 = 9.1875;            // diameter of the inner hole
	double d2 = 1.0088;            // width of the ring hole
	double L = 15;                 // wavelength of the soundwave
	double H = 12.8219;            // viewing distance
	double hole_radius = 59.0000;  // output wave radius
						   /*
						   there is a flat sound absorbing surface
						   inside that surface there is a circular hole with diameter d1
						   also there is a ring - shaped hole with inner radius R - d2 / 2 and outer radius R + d2 / 2
						   the inner hole is positioned in the center of the ring hole
						   also there is a point sound source
						   */

	double goodness;
	double contrast;
	double size_allignment;
	double size_constancy;
	double output_power;
	double output_dissipation;
	bool has_a_zero_point;

	double minimal_intensity = INT32_MAX;
	double maximal_intensity = 0;
	double wave_r1 = 0; // output wave radius on double distance

	Shape() {
		hole_radius = 50;
	}
	void print_data() {
		calc_goodness();

		cout << endl;
		cout << "------------------------------------------------------------------------" << endl;
		cout << "Viewing height:                           " << H << endl;
		cout << endl;
		cout << "Inner hole diameter:                      " << d1 << endl;
		cout << "Ring hole radius:                         " << R << endl;
		cout << "Ring hole width:                          " << d2 << endl;
		cout << "Output hole diameter:                     " << hole_radius * 2 << endl;
		cout << "Wavelength:                               " << L << endl;
		cout << endl;
		cout << "Minimal sound intensity:                  " << minimal_intensity << endl;
		cout << "Maximal sound intensity:                  " << maximal_intensity << endl;
		cout << "Output wave diameter:                     " << get_radius(H, 1) * 2 << endl;
		cout << "Output wave diameter on double distance:  " << wave_r1 * 2 << endl;
		cout << endl;
		cout << "Relative goodness:                        " << goodness << endl;
		cout << "Contrast:                                 " << contrast << endl;
		cout << "Output wave size allignment:              " << size_allignment << endl;
		cout << "Output wave size constancy:               " << size_constancy << endl;
		cout << "Wave after wall power:                    " << output_power << endl;
		cout << "Wave after wall dissipaton:               " << output_dissipation << endl;
		cout << "Has a minimum:                            " << has_a_zero_point << endl;
		cout << "------------------------------------------------------------------------" << endl;
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
	void optimize() {
		cout << endl << "Optimisation started" << endl << endl;

		cout << "Start data:" << endl;
		print_data();

		for (double optimisation_coefficient = 1; optimisation_coefficient > 0.05; optimisation_coefficient -= 0.05) {
			find_local_maximum(25, optimisation_coefficient, 0, 0);
		}

		cout << "End data:" << endl;
		print_data();

		cout << endl << "Optimisation ended" << endl << endl;
	}
	bool check_validity() {
		if (H < L / 2) {
			return 0;
		}
		if (R - d2 / 2 <= d1 / 2) {
			return 0;
		}

		return 1;
	}
	double get_radius(double h, bool globalize) {
		double out;

		double minimal = INT32_MAX;
		double maximal = 0;
		bool zero_point = 0;
		double old = getintensity(XSize / 2, YSize / 2, h);
		for (ush new_r = 1; new_r < XSize / 2; new_r++) {
			double intensity = getintensity(XSize / 2, YSize / 2 + new_r, h);

			if (intensity > old) {
				if (intensity < minimal) {
					out = new_r;
					zero_point = 1;
					minimal = old;
				}
			}
			if (intensity > maximal) {
				maximal = intensity;
			}

			old = intensity;
		}

		if (globalize) {
			minimal_intensity = minimal;
			maximal_intensity = maximal;
			has_a_zero_point = zero_point;
		}

		return out;
	}
	double show_slice() {
		minimal_intensity = INT32_MAX;
		maximal_intensity = 0;

		byte color[3];
		double intensity;
		for (ush xcor = 0; xcor < XSize; xcor++) {
			for (ush ycor = 0; ycor < YSize; ycor++) {
				intensity = getintensity(xcor, ycor, H);

				if (intensity > maximal_intensity) {
					maximal_intensity = intensity;
				}
				if (intensity < minimal_intensity) {
					minimal_intensity = intensity;
				}

				color[0] = 127 * (1 + sin(pow(intensity, 0.15) * pi * 0.5 + pi * 0));
				color[1] = 127 * (1 + sin(pow(intensity, 0.15) * pi * 1 + pi * 0));
				color[2] = 117 * (1 + sin(pow(intensity, 0.15) * pi * 0.5 + pi * 0.5));

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
		byte color[3];
		double intensity;
		for (ush h = 1; h < XSize + 1; h++) {
			for (ush ycor = 0; ycor <= YSize / 2; ycor++) {
				intensity = getintensity(XSize / 2, ycor, h);

				color[0] = 127 * (1 + sin(pow(intensity, 0.15) * pi * 0.5 + pi * 0));
				color[1] = 127 * (1 + sin(pow(intensity, 0.15) * pi * 1 + pi * 0));
				color[2] = 117 * (1 + sin(pow(intensity, 0.15) * pi * 0.5 + pi * 0.5));

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
		get_radius(H, 1);

		for (ush r = 0; r < hole_radius; r++) {
			intensity_profile[r] = getintensityonwall(XSize / 2, YSize / 2 + r);
		}

		byte color[3];
		double intensity;
		for (ush d = 1; d < XSize + 1; d++) {
			for (ush ycor = 0; ycor <= YSize / 2; ycor++) {
				intensity = getintensityafterwall(XSize / 2, ycor, d, intensity_profile);

				color[0] = 127 * (1 + sin(pow(intensity, 0.15) * pi * 0.5 + pi * 0));
				color[1] = 127 * (1 + sin(pow(intensity, 0.15) * pi * 1 + pi * 0));
				color[2] = 117 * (1 + sin(pow(intensity, 0.15) * pi * 0.5 + pi * 0.5));

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
	cout << "shound get a direct soundwave. The distance from the hole to the barier is H." << endl;
	cout << endl;
	cout << "List of commands:" << endl;
	cout << "    help" << endl;
	cout << "        shows help" << endl << endl;
	cout << "    exit" << endl;
	cout << "        breaks the program" << endl << endl;
	cout << "    show data" << endl;
	cout << "        shows shape data" << endl << endl;
	cout << "    show slice" << endl;
	cout << "        shows wave slice on chosen H" << endl << endl;
	cout << "    show profile" << endl;
	cout << "        shows wave profile" << endl << endl;
	cout << "    show after wall" << endl;
	cout << "        shows wave profile after exiting the wall" << endl << endl;
	cout << "    optimize" << endl;
	cout << "        automatically optimizes the shape" << endl << endl;
	cout << "    H" << endl;
	cout << "        changes value of H" << endl << endl;
	cout << "    R" << endl;
	cout << "        changes value of R" << endl << endl;
	cout << "    r" << endl;
	cout << "        changes value of r" << endl << endl;
	cout << "    d1" << endl;
	cout << "        changes value of d1" << endl << endl;
	cout << "    d2" << endl;
	cout << "        changes value of d2" << endl << endl;
	cout << "    L" << endl;
	cout << "        changes value of L" << endl;
	cout << endl;
	cout << "Red color on wave visualisation means, that the amplitude of the wave is high" << endl;
	cout << "Green color on wave visualisation means, that the amplitude of the wave is medium" << endl;
	cout << "Red color on wave visualisation means, that the amplitude of the wave is low" << endl;
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
	else if (command == "show data") {
		main_shape.print_data();
	}
	else if (command == "show slice") {
		main_shape.show_slice();
		main_shape.print_data();
	}
	else if (command == "show profile") {
		main_shape.show_profile();
		main_shape.print_data();
	}
	else if (command == "show after wall") {
		main_shape.show_after_wall();
	}
	else if (command == "local maximum") {
		ush number_of_iterations;
		double optimisation_value;
		cout << "Enter number of iterations: ";
		cin >> number_of_iterations;
		cout << "Enter the optimisation value: ";
		cin >> optimisation_value;
		main_shape.find_local_maximum(number_of_iterations, optimisation_value, 1, 0.3);
	}
	else if (command == "optimize") {
		main_shape.optimize();
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