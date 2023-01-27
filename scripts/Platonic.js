//	Platonic.js
//	Requires animate.js, BSP.js, GonND.js, Gon3D.js, SphericalTriangles.js, Colours.js, GraphicsUtils.hs, Utils.js included first

w=0;
h=0;
wc=0;
hc=0;
bg=null;
F = 0;

resized=false;
dpr = 1;

//	Polyhedra

pair=-1;
npoly=2;
ph=new Array(npoly);
rad=0;

rotation=Utils.multiDim([npoly,3,3],0);

gons=null;
sb=null;
sf=null;

//	Rendering

nshades=3;
palette=null;
lum=null;
xp=null;
yp=null;

customAnimate.initApp1 =
function ()
{
dpr=Math.min(Math.floor(window.devicePixelRatio || 1),2);
if (dpr > 1 && !resized)
	{
	canvas.style.width = canvas.width+"px";
	canvas.style.height = canvas.height+"px";
	canvas.width *= dpr;
	canvas.height *= dpr;
	resized = true;
	};
	
w=canvas.width;
h=canvas.height;
ctx=canvas.getContext('2d');
wc=w/2;
hc=h/2;

if (bg==null) bg=ctx.getImageData(0,0,w,h);
ctx.lineWidth=dpr*0.75;

rad=Math.min(wc,hc)-4;

palette=GraphicsUtils.randomPalette(npoly+1,nshades,1);

palette[0][0]+="aa";
palette[0][1]+="aa";

//	Pick types of polyhedron (0=tetrahedron ... 4=icosahedron)

var indexList=[[0,0],[1,1],[1,2],[2,2],[2,1],[3,3],[3,4],[4,4],[4,3]];
var ni=indexList.length;
if (pair<0) pair=Math.floor(ni*Math.random());
else pair=(pair+1)%ni;
var index=indexList[pair];

var maxgons=0;
for (var ip=0;ip<npoly;ip++)
	{
	var phA=ph[ip]=Polyhedron.regularPolyhedron(index[ip]);
	maxgons+=4*phA.NF;

	//	Random rotation matrix
	
	Utils.setRotation(rotation[ip],
		0.02*(2.0+Math.random()),
		Utils.randomUnitVec(false),
		1.0);
	};

gons=new Array(maxgons);	//	Merged list of polygons, both unbroken
							//	original faces, and inside segments.

var p=ph[0];
var nf=p.NF, fv=p.face_vertices;
var vertices=p.vertices;

//	Put all faces of first polyhedron into Gon3D list

for (var jf=0;jf<nf;jf++)
	{
	var fvj=fv[jf];
	var ns=fvj.length;
	var vt=new Array(ns);
	for (var i=0;i<ns;i++) vt[i]=vertices[fvj[i]];
	gons[jf]=new Gon3D(vt,ns);
	gons[jf].surfaceID=jf;
	};
		
sb=new Array(nf);			//	Dot product of light source reflection
							//	with viewing vector
sf=new Array(nf);			//	True for front faces

//	Rendering

lum=Utils.randomUnitVec(false);
xp=Utils.multiDim([maxgons,50],0);
yp=Utils.multiDim([maxgons,50],0);

F = 0;
};

customAnimate.draw1 =
function (t)
{
ctx.putImageData(bg,0,0);

for (var ip=0;ip<npoly;ip++) ph[ip].transform(rotation[ip],1.0);

//	Get face intersections
		
var igons=Polyhedron.faceIntersections(ph,
	Polyhedron.POLYHEDRON_0 | Polyhedron.ALL_FACES | Polyhedron.INTERIOR_PIECE);

var ngons=0;
for (var ip=0;ip<1;ip++)				//	For each polyhedron
	{
	var phA=ph[ip];
	var nf=phA.NF;
	var normals=phA.faces;
	
	for (var jf=0;jf<nf;jf++)			//	For each face
		{
		//	Norm
		
		var norm=normals[jf];
		sf[jf]=norm[2]>0.0;
		
		//	Brightness calculations

		var ndl=0.0;
		for (var i=0;i<3;i++) ndl+=norm[i]*lum[i];
		var brightness=lum[2]-2.0*ndl*norm[2];
		sb[jf]=brightness;
		};
	ngons+=nf;
	};
	
//	Add interior pieces to merged gons list

var np=igons.length;
for (var i=0;i<np;i++) gons[ngons++]=igons[i];
	
//	Put all gons into x and y coordinate lists

for (var i=0;i<ngons;i++)
	{
	var gnd=gons[i];
	var gi=gnd.vertices;
	var gn=gnd.nvert;
	var xpi=xp[i], ypi=yp[i];
	for (var j=0;j<=gn;j++)
		{
		var gk=gi[j%gn];
		xpi[j]=(wc+rad*gk[0]);
		ypi[j]=(hc+rad*gk[1]);
		};
	};
	
//	Render
//	------

//	pass 0: Front-facing full faces
//	pass 1: Borders for back-facing full faces
//	pass 2: Back-facing interior segments
//	pass 3: Front-facing interior segments
//	pass 4: Borders for front-facing segments
//	pass 5:	Borders for front-facing full faces

for (var pass=0;pass<6;pass++)
	{
	for (var i=0;i<ngons;i++)
		{
		var gnd=gons[i];
		var gn=gnd.nvert;
		var fi=gnd.surfaceID;
		var brightness=sb[fi];
		var front=sf[fi], back=!front;
		var full=(gnd.flags==0), seg=!full;
		var gc=full ? 0 : (front ? 1 : 2);
		
		if (pass==0 && front && full ||
			pass==2 && back && seg ||
			pass==3 && front && seg)
			{
			ctx.fillStyle = palette[gc][back ? 0: (brightness<0.1?0:(brightness>0.6?2:1))];
			GraphicsUtils.polygon(ctx,xp[i],yp[i],gn,true,false);
			}
		else if (
			pass==1 && back && full ||
			pass==4 && front && seg ||
			pass==5 && front && full)
			{
			ctx.strokeStyle=palette[gc][front ? 0 : 2];
			GraphicsUtils.polygon(ctx,xp[i],yp[i],gn+1,false,true);
			};
		};
	};

F++;
};
