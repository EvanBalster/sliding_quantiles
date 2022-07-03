#include <iostream>
#include <iomanip>
#include <string>
#include <array>
#include <deque>

#include <quern/histogram_tracked.hpp>


using namespace quern::literals;


quern::quantile_fraction<> p_quantiles[] =
{
#if 0
	{1,2}
#else
	// Extrema
	1/100_quo, 5/100_quo, 10/100_quo,

	// Median & quartiles
	1/4_quo,
	1/2_quo,
		2/4_quo,
	3/4_quo,

	// Extrema
	90/100_quo, 95/100_quo, 99/100_quo,
#endif
};


using Histogram32 = quern::histogram<float>;


struct QuantileTester :
	public quern::histogram_tracked<Histogram32>
{
public:
	QuantileTester() :
		histogram_tracked(quern::binning_params<float>{0.f, 32.f, 32})
	{
		histogram_tracked::add_quantiles(p_quantiles);
	}

	~QuantileTester()
	{
	}

#if 1
	void print()
	{
		auto &hist = histogram();
		std::cout << "\tHistogram:  population " << population() << std::endl;
		for (auto i = hist.begin(), e = hist.end(); i < e; ++i)
		{
			for (auto &q : quantiles())
				if (q.index_range.is_range() && q.index_range.upper == i.index())
					std::cout << "\t\t\t<-" << q.quantile.num << "/" << q.quantile.den
						<< "  (" << q.index_range.lower << "," << q.index_range.upper << ")" << std::endl;

			if (*i == 0) continue;

			// Histogram bar
			std::cout << "\t" << std::setw(4) << i.index() << ":" << std::setw(5) << *i << " ";
			size_t nibs = 100 - (100*(population()-*i) / population());
			for (size_t i = 0; i < nibs; ++i)
			{
				std::cout << '=';
			}

			// Value quantiles at this position
			unsigned quantiles_here = 0;
			for (auto &q : quantiles())
				if (q.index_range.is_value() && q.index_range.lower == i.index())
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
					<< "): population " << population() << std::endl;
			printedHeading = true;
		};

		{
			auto &hist = histogram();

			// Population verification
			{
				auto expect = hist.calc_population();
				if (population() != expect)
				{
					printHeading();
					std::cout << "\t\tPopulation inconsistent: actual total is " << expect
						<< " but cached value is " << population() << std::endl;
				}
			}

			// Correct samples_lower
			for (auto &q : quantiles())
			{
				size_t count = 0;
				for (size_t i = 0, e = q.index_range.upper; i < e; ++i)
					count += hist.count_at(i);

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
				auto expected = find_quantile_indexes(hist, q.quantile);

				bool consistent = true;
				if (expected.lower != q.index_range.lower) consistent = false;
				if (expected.upper != q.index_range.upper) consistent = false;

				if (!consistent)
				{
					printHeading();
					std::cout << "\t\tBad quantile " << q.quantile.num << "/" << q.quantile.den
						<< " .. location is " << q.index_range.lower << ":" << q.index_range.upper
						<< " but histogram evaluates to " << expected.lower << ":" << expected.upper
						<< " ... last adjust " << q.last_adjust
						<< std::endl;
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
					<< " " << std::setw(3) << q.index_range.lower << ":"
					<< std::left << std::setw(3) << q.index_range.upper << std::right
					<< " samples_lower = " << q.samples_lower
					<< " ... last adjust " << q.last_adjust
					<< std::endl;
			}

			std::cout << "\t**********" << std::endl;
		}
#endif
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
				for (size_t maxroll = 0; maxroll < n; maxroll += 3+(maxroll&1)) x += rand() % (4+(maxroll&1));
				test.insert(x);
				test.consistencyCheck("random insertion");
			}
			test.print();
		}
	}

	for (size_t pop = 10; pop < 10000; pop = pop*3 + pop/2)
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
