/*

Gon3D.js
==========

Class representing a planar polygon in 3 dimensions.

This specialises the class GonND, for n-dimensional polygons.

*/

//	Constructor
//	===========

function Gon3D()
{
var narg = arguments.length;

//	Use appropriate GonND constructor

if (narg==2) GonND.call(this,3,arguments[0],arguments[1]);
else if (narg==3) GonND.call(this,3,arguments[0],arguments[1],arguments[2]);
else if (narg==4) GonND.call(this,3,arguments[0],arguments[1],arguments[2],arguments[3]);
else
	{
	throw "Bad constructor for GonND"
	};
}
Gon3D.prototype = GonND.prototype;


//	newGon()
//	========

//	newGon() returns a new GonND by copying vertex pointers;
//	subclasses override this, to return a GonND with the full specificity of
//	the subclass.

Gon3D.prototype.newGon =
function (v, s, n)
{
return new Gon3D(v, s, n);
};


//	plane()
//	=======
//
//	Create a 4-dimensional array which describes the plane the polygon lies in.
//
//	plane[0,1,2]	is a unit normal to the plane.  Direction is chosen so that
//					looking back (i.e. unit normal towards viewer) polygon is traversed
//					in a counter-clockwise direction.
//	plane[3]		is the dot product of that unit normal with any point on the
//					plane, i.e. the distance of the plane from the origin along
//					the orthogonal direction.

Gon3D.prototype.plane =
function ()
{
var p=[0,0,0,0];

//	Scan for a pair of non-collinear edges

var ssq=0.0, dp=0.0, pi;
for (var bv=0;bv<this.nvert;bv++)
	{
	var v0=this.vertices[bv], v1=this.vertices[(bv+1)%this.nvert], v2=this.vertices[(bv+2)%this.nvert];
	ssq=dp=0.0;
	for (var i=0;i<3;i++)
		{
		var j=(i+1)%3;
		var k=(i+2)%3;
		p[i]=pi=(v1[j]-v0[j])*(v2[k]-v1[k])-(v1[k]-v0[k])*(v2[j]-v1[j]);
		ssq+=pi*pi;
		dp+=pi*v0[i];
		};
	if (ssq>1e-8) break;		//	Edges not collinear
	};
ssq=Math.sqrt(ssq);
for (var i=0;i<3;i++) p[i]/=ssq;
p[3]=dp/ssq;
return p;
}

//	markEdges()
//	===========

//	Split and/or mark those edges of the Gon3D that lie within one of the listed
//	objects, by setting the flag GonND.EDGE_IN_CONTENT for their edgeIDs entry.

Gon3D.prototype.markEdges =
function (svec, sContents, idToCheck, tolerance, pm)
{
var nc=0;
if (sContents==null || (nc=sContents.length)==0) return;

var ourPlane=this.plane();	//	Get plane for this Gon3D
var x=[0,0];
var flaggedID=idToCheck | GonND.EDGE_IN_CONTENT;

//	Process all objects in contents vector; these should be Gon3Ds

for (var ic=0;ic<nc;ic++)
	{
	var cobj=sContents[ic];
	if (!(cobj instanceof Gon3D)) continue;
	var cgon=cobj;
	
	//	Get information about the line segments where our plane intersects it
	
	var segs=cgon.splitSegs(ourPlane, tolerance, pm, false, null);
	var nsegs=segs.length;
	
	for (var iseg=0;iseg<nsegs;iseg++)
		{
		var p0=segs[iseg][0], p1=segs[iseg][1], segPts=[p0,p1];
		var segVec=[0,0,0];
		for (var i=0;i<3;i++) segVec[i]=p1[i]-p0[i];
		var segLen=Utils.normalise(segVec);
		var baseDP=Utils.dot(segVec,p0);
		
		//	For each split-created edge

		for (var ie=0;ie<this.nvert;ie++)
		if (this.edgeIDs[ie]==idToCheck)
			{
			//	Find its vertex parameters along segment
			
			var beforeSeg=0, bLow=0, inSeg=0, bHigh=0, afterSeg=0;
			for (var i=0;i<2;i++)
				{
				var xx=(Utils.dot(segVec,this.vertices[(ie+i)%this.nvert])-baseDP)/segLen;
				if (Math.abs(xx)<tolerance) {xx=0.0; bLow++;}
				else if (Math.abs(1.0-xx)<tolerance) {xx=1.0; bHigh++;}
				else if (xx<0) beforeSeg++;
				else if (xx>1.0) afterSeg++;
				else inSeg++;
				x[i]=xx;
				};
				
			//	Edge wholly in segment
			
			if (beforeSeg+afterSeg==0) this.edgeIDs[ie]=flaggedID;
			
			//	Edge wholly outside segment
			
			else if (bLow+beforeSeg==2 || bHigh+afterSeg==2) continue;
			
			//	Edge part in / part out, so we need to create
			//	new vertices and edges
			
			else
				{
				var nextra=beforeSeg+afterSeg;	//	Number of extra vertices
				var newIDs=new Array(this.nvert+nextra);
				var newVertices=new Array(this.nvert+nextra);
				
				//	Copy all existing vertices and edge IDs, and insert
				//	new ones along edge ie
				
				var newV=0;
				for (var i=0;i<this.nvert;i++)
					{
					newVertices[newV]=this.vertices[i];
					newIDs[newV++]=this.edgeIDs[i];
					
					if (i==ie)
						{
						var sp=0;			//	Count endpoints of segment we've used
						var inContent;
						if (inContent=(x[0]>=0.0 && x[0]<=1.0))
							{
							//	Edge is already in content from start, and
							//	we've passed one endpoint of segment
							
							newIDs[newV-1]=flaggedID;
							sp++;
							};
						
						for (var j=0;j<nextra;j++)
							{
							inContent=!inContent;
							newVertices[newV]=segPts[x[1]>x[0]?sp:1-sp];
							sp++;
							newIDs[newV++]=inContent?flaggedID:idToCheck;
							};
						};
					};
				
				ie+=nextra;
				this.nvert+=nextra;
				this.edgeIDs=newIDs;
				this.vertices=newVertices;
				};
			};
		};
	};
}

