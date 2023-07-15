#include <fstream>

#include "Image.h"

using namespace std;


int main()
{
	Image img(1280,640);
	
	img.penColor(25,25,25);
	img.text(210,100, "single header library", IBM_EGA_8x14, 5);
	img.penColor(12,12,12);
	img.text(170,150, "MS-DOS", IBM_EGA_8x14, 20);
	
	
	img.penColor(255,255,255);
	img.text(510,250, "SLAP", Portfolio_6x8, 12);
	img.text(300,400, "Simple Linear Algebra Package", ISO_font, 3);
	
	img.penColor(100,100,100);
	int N = 20;
	int Delta_x = img.width()/N;
	int Delta_y = img.height()/N;
	for(int i=0; i<=N; i++){
		int x = i * Delta_x;
		int y = i * Delta_y;
//		cout << x << " " << y << endl;
		img.line(0, y-1, x, img.height()-1);
	}
	img.penColor(255,255,255);
	img.penWidth(2);
	img.rect(430,200, 850,380);
	img.penWidth(3);
	img.rect(1,1,img.width()-3,img.height()-2);
	
	img.save_bmp("cover.bmp");
}
