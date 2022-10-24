#include "pbPlots.hpp"
#include "supportLib.hpp"

using namespace std;

double* range(int in, int out) {
	double* res = new double[out + 1];
	for (int i = in; i <= out; i++) {
		res[i] = i;
	}
	return res;
}

bool plot(double* x, double* y, int in, int out, string fileName) {
	bool success;
	StringReference* errorMessage = CreateStringReferenceLengthValue(0, L' ');
	RGBABitmapImageReference* imageReference = CreateRGBABitmapImageReference();

	vector<double> xs(x + in, x + out + 1);
	vector<double> ys(y + in, y + out + 1);

	double xMin = x[in], xMax = x[in], yMin = y[in], yMax = y[in];
	for (int i = in + 1; i <= out; i++) {
		if (x[i] < xMin) xMin = x[i];
		if (x[i] > xMax) xMax = x[i];
		if (y[i] < yMin) yMin = y[i];
		if (y[i] > yMax) yMax = y[i];
	}
	double xRange = xMax - xMin, yRange = yMax - yMin;
	const double RELATIVE_MARGIN = 0.03;

	ScatterPlotSeries* series = GetDefaultScatterPlotSeriesSettings();
	series->xs = &xs;
	series->ys = &ys;
	series->linearInterpolation = false;
	series->pointType = toVector(L"crosses");
	series->color = GetGray(0.9);
	//series->color = CreateRGBColor(1, 0, 1);

	ScatterPlotSettings* settings = GetDefaultScatterPlotSettings();
	settings->width = 900;
	settings->height = 600;
	settings->autoBoundaries = false;
	settings->xMin = xMin - xRange * RELATIVE_MARGIN;
	settings->xMax = xMax + xRange * RELATIVE_MARGIN;
	settings->yMin = yMin - yRange * RELATIVE_MARGIN;
	settings->yMax = yMax + yRange * RELATIVE_MARGIN;
	settings->autoPadding = true;
	settings->title = toVector(L"");
	settings->xLabel = toVector(L"");
	settings->yLabel = toVector(L"");

	settings->scatterPlotSeries->push_back(series);

	success = DrawScatterPlotFromSettings(imageReference, settings, errorMessage);

	if (success) {
		vector<double>* pngdata = ConvertToPNG(imageReference->image);
		WriteToFile(pngdata, fileName);
		DeleteImage(imageReference->image);
	}
	else {
		cerr << "Error: ";
		for (int i = 0; i < errorMessage->string->size(); i++) {
			wcerr << errorMessage->string->at(i);
		}
		cerr << endl;
	}

	FreeAllocations();

	return success ? 0 : 1;
}