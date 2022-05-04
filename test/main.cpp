#include <iostream>
#include <iomanip>
#include <string>
#include <array>
#include <deque>

#include <histogram_quantiles.hpp>


edb::quantile_position<> p_quantiles[] = {

#if 0
	{{1,2}}
#else
	// Extrema
	{{1,100}}, {{5,100}}, {{10,100}},

	// Median & quartiles
	{{1,4}},
	{{1,2}},
		{{2,4}},
	{{3,4}},

	// Extrema
	{{90,100}}, {{95,100}}, {{99,100}},
#endif
};


using Histogram32 = edb::histogram_basic<std::array<size_t, 32>>;

struct QuantileTester :
	public edb::quantile_tracker<Histogram32>
{

	QuantileTester() :
		quantile_tracker(std::vector<edb::quantile_position<>>
			(p_quantiles, p_quantiles+ sizeof(p_quantiles)/sizeof(p_quantiles[0]) ))
	{
	}

	~QuantileTester()
	{
	}

	void print()
	{
		auto &hist = histogram();
		std::cout << "\tHistogram:  population " << hist.population() << std::endl;
		for (size_t i = 0; i < hist.size(); ++i)
		{
			for (auto &q : quantiles())
				if (q.range.is_range() && q.range.upper == i)
					std::cout << "\t\t\t<-" << q.quantile.num << "/" << q.quantile.den
						<< "  (" << q.range.lower << "," << q.range.upper << ")" << std::endl;

			if (hist[i] == 0) continue;

			// Histogram bar
			std::cout << "\t" << std::setw(4) << i << ":" << std::setw(5) << hist[i] << " ";
			size_t nibs = 100 - (100*(hist.population()-hist[i]) / hist.population());
			for (size_t i = 0; i < nibs; ++i)
			{
				std::cout << '=';
			}

			// Value quantiles at this position
			unsigned quantiles_here = 0;
			for (auto &q : quantiles())
				if (q.range.is_value() && q.range.lower == i)
					std::cout << (quantiles_here++ ? ", " : " <- ") << q.quantile.num << "/" << q.quantile.den;

			std::cout << std::endl;
		}
	}

	void consistencyCheck(std::string context)
	{
		bool printedHeading = false;
		auto printHeading = [&]()
		{
			if (!printedHeading)
				std::cout << "\tConsistency Checks (" << context
					<< "): population " << histogram().population() << std::endl;
			printedHeading = true;
		};

		{
			auto &hist = histogram();

			// Correct samples_lower
			for (auto &q : quantiles())
			{
				size_t count = 0;
				for (size_t i = 0, e = q.range.upper; i < e; ++i)
					count += hist[i];

				if (count != q.samples_lower)
				{
					printHeading();
					std::cout << "\t\tInconsistency at " << q.quantile.num << "/" << q.quantile.den
						<< " samples_lower is " << q.samples_lower
						<< " but should be " << count << std::endl;
				}
			}


			// Correct quantile values
			for (auto &q : quantiles())
			{
				auto expected = hist.find_quantile(q.quantile.num, q.quantile.den);

				bool consistent = true;
				if (expected.lower != q.range.lower) consistent = false;
				if (expected.upper != q.range.upper) consistent = false;

				if (!consistent)
				{
					printHeading();
					std::cout << "\t\tBad quantile " << q.quantile.num << "/" << q.quantile.den
						<< " .. location is " << q.range.lower << ":" << q.range.upper
						<< " but histogram evaluates to " << expected.lower << ":" << expected.upper << std::endl;
				}
			}
		}

		if (printedHeading)
		{
			print();

			std::cout << "\tQuantile data:" << std::endl;
			for (auto &q : quantiles())
			{
				std::cout
					<< "\t\t" << std::setw(3) << q.quantile.num
					<< "/" << std::left << std::setw(3) << q.quantile.den << std::right
					<< " " << std::setw(3) << q.range.lower << ":"
					<< std::left << std::setw(3) << q.range.upper << std::right
					<< " samples_lower = " << q.samples_lower
					<< " ... last adjust " << q.last_adjust
					<< std::endl;
			}

			std::cout << "\t**********" << std::endl;
		}
	}
};


int main(int argc, char **argv)
{
	std::srand(clock());

	for (size_t n = 2; n < 20; n += (1+n/4))
	{
		{
			std::cout << "TEST: Rectangular up to " << n << std::endl;

			{
				QuantileTester test;
				for (size_t i = 0; i < n; ++i)
				{
					test.insert(i);
					test.consistencyCheck("insertion");
				}
				//test.print();
			}
			std::cout << std::endl;
		}

		{
			std::cout << "TEST: Rectangular down to " << n << std::endl;
			{
				QuantileTester test;
				for (size_t i = n; i--;)
				{
					test.insert(i);
					test.consistencyCheck("insertion");
				}
				//test.print();
			}
			std::cout << std::endl;
		}
	}

	for (size_t n = 2; n <= 32; n *= 2)
	{
		std::cout << "TEST: 1000 random insertions over range " << n << std::endl;

		{
			QuantileTester test;
			for (size_t i = 0; i < 1000; ++i)
			{
				size_t x = 0;
				for (size_t bit = 0; bit < n; bit += 2) x += rand() % 4;
				test.insert(x);
				test.consistencyCheck("random insertion");
			}
			test.print();
		}
	}

	for (size_t pop = 10; pop < 1000; pop = pop*3 + pop/2)
	{
		std::cout << "TEST: rolling insertions, population " << pop << std::endl;

		{
			std::deque<size_t> log;

			QuantileTester test;

			for (size_t i = 0; i < pop; ++i)
			{
				size_t x = size_t(rand()) & 31;
				test.insert(x);
				log.push_back(x);
				test.consistencyCheck("rolling insertion, pre-fill");
			}

			test.print();

			for (size_t i = 0; i < 10000; ++i)
			{
				size_t x = size_t(rand()) & 31;
				//std::cout << "\t\treplace " << log.front() << " with " << x << std::endl;
				test.replace(x, log.front()); log.pop_front();
				log.push_back(x);
				test.consistencyCheck("rolling insertion, run");
			}

			test.print();

			while (log.size())
			{
				test.remove(log.front()); log.pop_front();
				test.consistencyCheck("rolling insertion, empty out");
			}

			test.print();
		}
	}

	std::cout << "Complete.  Press ENTER to close." << std::endl;
	std::cin.ignore(255, '\n');
}
