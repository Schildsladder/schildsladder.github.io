//	GraphicsUtils.js
//	Utility graphics routines

GraphicsUtils = {};


//	Approximate Gaussian integral for drawing pixellated Gaussian disks
//	===================================================================
//
//	This function returns:
//
//	0																x <= -1
//	a polynomial approximation of the integral of a Gaussian		-1 < x < 1
//	1																x >= 1

GraphicsUtils.approximateGaussianIntegral =
function (x)
{
if (x <= -1.0) return 0.0;
if (x >= 1.0) return 1.0;
var x2=x*x;
return 0.5+(945.0 - x2*(840.0 - x2*(504.0 - x2*(180.0 - 35.0*x2))))*x/928.0;
};

//	Return a set of pixellated brightnesses for a 1-dimensional Gaussian bump of specified centre and radius
//	========================================================================================================

GraphicsUtils.pixellatedGaussian =
function (x0, r, res, result)
{
var lo=Math.floor(x0-r), hi=Math.floor(x0+r), n=(hi-lo+1);
var inc = 1.0/res;

var prev=0, max=GraphicsUtils.approximateGaussianIntegral(0.5*inc/r)-GraphicsUtils.approximateGaussianIntegral(-0.5*inc/r);
var rp=0;
for (var i=lo;i<=hi;i++)
for (var ii=1;ii<=res;ii++)
	{
	var x=(i+ii*inc-x0)/r;
	var g=GraphicsUtils.approximateGaussianIntegral(x);
	result[rp++]=(g-prev)/max;
	prev=g;
	};
return n;
}

//	Draw a Gaussian disk with specified centre, hue and saturation into an array of canvas image data
//	=================================================================================================

GraphicsUtils.a1 = null;
GraphicsUtils.a2 = null;
GraphicsUtils.rgb = new Array(3);

GraphicsUtils.gaussianDisk =
function (imgData, w, h, x0, y0, r, hue, sat, bri, res)
{
var rgb=GraphicsUtils.rgb, a1, a2;
var len=Math.ceil((2*r+2)*res);
var res2=res*res;

if (GraphicsUtils.a1==null || GraphicsUtils.a1.length<len) GraphicsUtils.a1=new Array(len);
if (GraphicsUtils.a2==null || GraphicsUtils.a2.length<len) GraphicsUtils.a2=new Array(len);

var nx=GraphicsUtils.pixellatedGaussian(x0, r, res, a1=GraphicsUtils.a1);
var ny=GraphicsUtils.pixellatedGaussian(y0, r, res, a2=GraphicsUtils.a2);

var loX=Math.floor(x0-r), loY=Math.floor(y0-r);
for (var ix=0;ix<nx;ix++)
	{
	var xx=loX+ix;
	if (xx>=0 && xx<w)
		{
		for (var iy=0;iy<ny;iy++)
			{
			var yy=loY+iy;
			if (yy>=0 && yy<h)
				{
				var z=0;
				for (var ii=0;ii<res;ii++)
					{
					var ax = a1[ix*res+ii];
					for (var jj=0;jj<res;jj++) z += ax * a2[iy*res+jj];
					};
				
				Colours.HSBtoRGB(hue,sat,(z/res2)*bri,rgb);
				var p = 4*(w*yy+xx);
				imgData[p] += 255*rgb[0];
				imgData[p+1] += 255*rgb[1];
				imgData[p+2] += 255*rgb[2];
				};
			};
		};
	};
};

//	Draw an ellipse
//	===============

GraphicsUtils.drawEllipse =
function (context, cx, cy, rx, ry, fill, stroke)
{
context.save();
context.beginPath();
context.translate(cx, cy);
context.scale(rx, ry);
context.arc(0, 0, 1, 0, 2 * Math.PI, false);
context.restore();
if (fill) context.fill();
if (stroke) context.stroke();
}

//	Status panel attached to an element
//	===================================

GraphicsUtils.statusPanel=null;

GraphicsUtils.showStatus =
function (el, txt)
{
var p=GraphicsUtils.associatePanel(el);
p.innerHTML = txt;
p.style.display="block";
}

GraphicsUtils.hideStatus =
function (el)
{
var p=GraphicsUtils.associatePanel(el);
p.innerHTML = "";
p.style.display="none";
}


GraphicsUtils.associatePanel =
function (el)
{
var newPanel = (GraphicsUtils.statusPanel == null);
if (newPanel) GraphicsUtils.statusPanel = document.createElement("p");
var st=GraphicsUtils.statusPanel.style;
if (newPanel)
	{
	st.position="absolute";
	st.zIndex=1;
	st.backgroundColor="#ffffff";
	st.color="#000000";
	st.padding="5px";
	st.margin="0px";
	};
var elRect = el.getBoundingClientRect();
var x = (window.pageXOffset !== undefined)
	? window.pageXOffset
	: (document.documentElement || document.body.parentNode || document.body).scrollLeft;
var y = (window.pageYOffset !== undefined)
	? window.pageYOffset
	: (document.documentElement || document.body.parentNode || document.body).scrollTop;
st.top=Math.round(elRect.top+y)+"px";
st.left=Math.round(elRect.left+x)+"px";

if (newPanel) document.body.appendChild(GraphicsUtils.statusPanel);
return GraphicsUtils.statusPanel;
}


//	A crude seeded random number generator
//	See:  https://stackoverflow.com/questions/521295/seeding-the-random-number-generator-in-javascript

GraphicsUtils.randomSeed = 1;

GraphicsUtils.random =
function ()
{
var x = Math.sin(GraphicsUtils.randomSeed++) * 10000;
return x - Math.floor(x);
}

//	Random multi-hue, multi-brightness palette
//	==========================================
//
//	tintCode is:
//
//	0	-	generate ncols completely independent colours
//	1	-	colours 2, 4, 6 ... are colours 1, 3, 5 ... filtered through colour 0
//	2	-	colours 1, 2, 3 ... are all filtered through colour 0

GraphicsUtils.randomPalette =
function (ncols, nshades, tintCode)
{
palette=new Array(ncols);
for (var i=0;i<ncols;i++) palette[i]=new Array(nshades);

//	HSB hues

var hue0=Math.random();
var hues=new Array(ncols);
for (var i=0;i<ncols;i++) hues[i]=(hue0+i/ncols)%1;

//	Generate RGB multiplicative filters to apply to light levels,

var filters=new Array(ncols), rgb=[0,0,0];
for (var icol=0;icol<ncols;icol++)
	{
	filters[icol]=[0,0,0];
	Colours.HSBtoRGB(hues[icol],1,1,rgb);
	var fil;
	for (var bgr=0;bgr<3;bgr++)
		{
		fil=rgb[2-bgr];
		if (icol>0)
			{
			if (tintCode==1)
				{
				if (icol%2==0) fil=filters[0][bgr]*filters[icol-1][bgr];
				}
			else if (tintCode==2) fil=(2+filters[0][bgr])/3*fil;
			};
		filters[icol][bgr]=fil;
		};

	//	Generate absolute RGB at various amounts of white light

	for (var ishade=0;ishade<nshades;ishade++)
		{
		var light=255-(nshades-1-ishade)*51;
		var rgbCol=0;
		for (var bgr=0;bgr<3;bgr++)
			{
			var col=Math.floor(filters[icol][bgr]*light);
			rgbCol|=(col << (bgr*8));
			};
		palette[icol][ishade]=Colours.RGBIntToStyle(rgbCol);
		};
	};

return palette;
}

//	polygon()
//	=========

GraphicsUtils.polygon =
function (ctx, xp, yp, np, fill, draw)
{
ctx.beginPath();
ctx.moveTo(xp[0],yp[0]);
for (var i=1;i<np;i++) ctx.lineTo(xp[i],yp[i]);
ctx.closePath();
if (fill) ctx.fill();
if (draw) ctx.stroke();
};


//	drawLine()
//	==========

GraphicsUtils.drawLine =
function (c,x1,y1,x2,y2)
{
c.beginPath();
c.moveTo(x1,y1);
c.lineTo(x2,y2);
c.stroke();
};

