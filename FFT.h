#include <memory>
#include <math.h>
#include <iterator>

static const double TWPI = 6.283185307179586;

class FFT
{
	int i, j, n, m, MAX, ISTEP;
	double Wpr, Wpi, Wr, Wi;
	double TMPR, TMPI, WTMP, THETA;
	std::unique_ptr<double[]> TMVL;
	std::unique_ptr<double[]> OUT;
	inline void calculate();
       inline void transform();
public:
	FFT(std::unique_ptr<double[]>in, int VALC);
	inline std::unique_ptr<double[]> get_out();
};
FFT::FFT(std::unique_ptr<double[]>in, int VALC)
{
	n = VALC * 2;
	m = VALC;
	std::unique_ptr<double[]>TMVL(new double[n]);
	for (i = 0; i < n; i += 2) {
		TMVL[i] = 0;
		TMVL[i + 1] = in[i / 2];//для рассчитывания поворотных множителей, оптимизация, берем значения точек только с четными индексами
	}
	initialize();
}
inline void FFT::initialize()
{
	MAX = 2;
	i = 1; j = 1;
	while (i < n) {
		if (j > i) {
			TMPR = TMVL[i];
			TMVL[i] = TMVL[j];
			TMVL[j] = TMPR;
			TMPR = TMVL[i + 1]; 
			TMVL[i + 1] = TMVL[j + 1];
			TMVL[j + 1] = TMPR;
		}
		i = i + 2; 
		while ((m >= 2) && (j > m)) {
			j = j - m;
			m = m >> 1;
		}
		j = j + m;
	}

	
	while (n > MAX) {
		THETA = -TWPI / MAX;
		Wpi = sin(THETA);
		WTMP = sin(THETA / 2);
		Wpr = WTMP * WTMP * 2;
		ISTEP = MAX * 2; Wr = 1; Wi = 0; m = 1;

		while (m < MAX) {
			i = m; m = m + 2; TMPR = Wr; TMPI = Wi;
			Wr = Wr - TMPR * Wpr - TMPI * Wpi;
			Wi = Wi + TMPR * Wpi - TMPR * Wpr;

			while (i < n) {
				j = i + MAX;
				TMPR = Wr * TMVL[j] - Wi * TMVL[j - 1];
				TMPI = Wi * TMVL[j] + Wr * TMVL[j - 1];

				TMVL[j] = TMVL[i] - TMPR; TMVL[j - 1] = TMVL[i - 1] - TMPI;//выразим елементы последовательностей
				TMVL[i] = TMVL[i] + TMPR; TMVL[i - 1] = TMVL[i - 1] + TMPI;
				i = i + ISTEP;
			}
		}

		MAX = ISTEP;
	}

}
inline void FFT::transform()
{
	std::unique_ptr<double[]> OUT(new double[m]);
	 for (i = 0; i < m; i++) {
	j = i * 2;
	OUT[i] = 2*sqrt(pow(TMVL[j],2) + pow(TMVL[j+1],2))/m;
  }
}
inline std::unique_ptr<double[]> FFT::get_out()
{
	std::cout << "in3:" << OUT[2] << "\n";//no debug
	return std::move(OUT);
}
