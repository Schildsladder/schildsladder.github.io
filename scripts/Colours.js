//	Colours.js
//	Utility colour routines

Colours = {};

//	Convert an HSB colour (with hue, sat, bri ranging from 0 to 1)
//	to an RGB triple (again scaled from 0 to 1)

Colours.HSBtoRGB =
function (hue,sat,bri,RGB)
{
var h6A = 6.0*hue;
var h6 = h6A - Math.floor(h6A/6);
var hi=Math.floor(h6);
var f=h6-hi;
var p=bri*(1.0-sat);
var q=bri*(1.0-f*sat);
var t=bri*(1.0-(1.0-f)*sat);
rgb = RGB || new Array(3);

if (hi==0)
	{
	rgb[0]=bri;
	rgb[1]=t;
	rgb[2]=p;
	}
else if (hi==1)
	{
	rgb[0]=q;
	rgb[1]=bri;
	rgb[2]=p;
	}
else if (hi==2)
	{
	rgb[0]=p;
	rgb[1]=bri;
	rgb[2]=t;
	}
else if (hi==3)
	{
	rgb[0]=p;
	rgb[1]=q;
	rgb[2]=bri;
	}
else if (hi==4)
	{
	rgb[0]=t;
	rgb[1]=p;
	rgb[2]=bri;
	}
else
	{
	rgb[0]=bri;
	rgb[1]=p;
	rgb[2]=q;
	};
	
return rgb;
};

//	Convert an HSB colour (with hue, sat, bri ranging from 0 to 1)
//	to a JavaScript colour style:  a hexadecimal string of the form "#RRGGBB"

Colours.HSBtoStyle =
function (hue,sat,bri)
{
return Colours.RGBtoStyle(Colours.HSBtoRGB(hue,sat,bri,null),255);
};

Colours.HSBtoRGBInt =
function(hue,sat,bri)
{
var rgb = Colours.HSBtoRGB(hue,sat,bri,null);
rgbInt = 0;
for (var i=0;i<3;i++) rgbInt = (rgbInt << 8) + Math.round(255*rgb[i]);
return rgbInt;
};

//	Convert an RGB colour triplet (with elements multiplied by f to get a range 0-255)
//	to a JavaScript colour style:  a hexadecimal string of the form "#RRGGBB"

Colours.RGBtoStyle =
function (rgb,f)
{
var res="#";
for (var i=0;i<3;i++)
	{
	var c=Math.round(rgb[i]*f);
	if (c==0) res+="00";
	else
		{
		if (c<16) res+="0";
		res+=(c).toString(16);
		};
	};
return res;
};

Colours.RGBAtoStyle =
function (rgba,f)
{
return "rgba("+Math.round(rgba[0]*f)+","+Math.round(rgba[1]*f)+","+Math.round(rgba[2]*f)+","+rgba[3]+")";
};

//	Convert an RGB integer 0xRRGGBB
//	to a JavaScript colour style:  a hexadecimal string of the form "#RRGGBB"

Colours.RGBIntToStyle =
function (rgbInt)
{
var res="";
for (var i=0;i<3;i++)
	{
	var c=rgbInt & 255;
	if (c==0) res="00"+res;
	else
		{
		res = (c).toString(16) + res;
		if (c<16) res="0"+res;
		};
	rgbInt = rgbInt >> 8;
	};
return "#"+res;
};

//	The  seven colours of the spectrum
//	==================================

Colours.spectrum=[
	[1.0, 0.0, 0.0],
	[1.0, 0.6, 0.0],
	[1.0, 1.0, 0.0],
	[0.0, 1.0, 0.0],
	[0.0, 0.0, 1.0],
//	[0.2, 0.4, 0.6],
	[0.2/0.6, 0.4/0.6, 1.0],
	[0.8, 0.6, 1.0]];
	
Colours.spectrumStyle =
function (i)
{
return Colours.RGBtoStyle(Colours.spectrum[i],255);
};

