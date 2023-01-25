	/*

SphericalTriangle.js
====================

Class representing a spherical triangle whose vertex angles are all rational
submultiples of PI.

*/

/*

Original Java class members:

//	The three angles at the vertices of the triangle are:
//
//		PI vA[i][0]/vA[i][1], i=0,1,2

public int vA[][];

//	The vertices (unit vectors)

public double vertices[][];

//	The cosines and sines of angular side lengths

private double sideCos[], sideSin[];

//	Triangles that own these vertices, and their owner vertex numbers

public SphericalTriangle vertexOwner[];
public int vertexNumber[];

//	Normals to the planes that define the edges of the triangle

public double edgeNormals[][];

//	The edges of the triangle can be linked to other triangles.

public SphericalTriangle links[];

//	If edge is linked, linkEdges gives the number of the link in the
//	neighbouring triangle.

public int linkEdges[];

//	The isometry that takes this triangle to the canonical position.

public double isometry[][];

//	Parity is +/-1 for copies/reflections of canonical triangle

public int parity;

//	Flags set during tiling which say we attached a new tile along the flagged
//	edge (i.e. bits 0, 1, 2 for edges 0, 1, 2)

public int tilingFlags;

//	Edge-link sequence that led to this triangle during tiling

public int nsteps, steps[];

//	User-defined ID

public int id;
*/

//	Constructors
//	------------


//customAnimate.draw1(t)

function SphericalTriangle()
{
var narg=arguments.length;

this.vA=null;
this.vertices=null;
this.sideCos=null;
this.sideSin=null;
this.vertexOwner=null;
this.vertexNumber=null;
this.edgeNormals=null;
this.links=null;
this.linkEdges=null;
this.isometry=null;
this.parity=0;
this.tilingFlags=0;
this.nsteps=0;
this.steps=null;
this.id=0;

if (narg==1)
	{

	//	Construct a spherical triangle in canonical position, with specified vertex angles.
	//
	//	Canonical position is:
	//
	//		vertex 0 is {0,0,1}
	//		vertex 1 is {sin(gamma), 0, cos(gamma)}, gamma to be determined
	//		vertex 2 is {cx,cy,cz}, all to be determined

	this.vA=arguments[0];
	this.links=new Array(3);
	this.linkEdges=new Array(3);

	//	Vertex angles

	var a=Math.PI*this.vA[0][0]/this.vA[0][1];
	var b=Math.PI*this.vA[1][0]/this.vA[1][1];
	var c=Math.PI*this.vA[2][0]/this.vA[2][1];

	//	Sines and cosines

	var ca=Math.cos(a), sa=Math.sin(a);
	var cb=Math.cos(b), sb=Math.sin(b);
	var cc=Math.cos(c), sc=Math.sin(c);

	//	Solve for canonical position

	var cz=(cb+ca*cc)/(sa*sc), omcz=Math.sqrt(Math.abs(1.0-cz*cz));
	var cx=ca*omcz, cy=Math.sqrt(Math.abs(1.0-cx*cx-cz*cz));
	var cosGamma=(Math.abs(ca)<1e-10 ? sa*cc/sb : (cz*sc-cb*sa)/(ca*sb));
	var sinGamma=sc*omcz/sb;

	var vt=[[0.0,0.0,1.0],[sinGamma,0.0,cosGamma],[cx,cy,cz]];
	this.vertices=vt;

	//	This triangle owns all its vertices

	this.vertexOwner=new Array(3);
	this.vertexNumber=new Array(3);
	for (var i=0;i<3;i++)
		{
		this.vertexOwner[i]=this;
		this.vertexNumber[i]=i;
		};

	//	Edge normals

	this.computeEdgeNormals();

	//	Isometry for canonical triangle is the identity

	this.isometry=Utils.multiDim([3,3],0);
	for (var i=0;i<3;i++) this.isometry[i][i]=1.0;
	this.parity=1;
	}
else
	{

	//	Construct a spherical triangle by reflection in one of the edges of another.
	//	Edge i is opposite vertex i, i=0,1,2.
	
	var neighbour = arguments[0];
	var edgeNumber = arguments[1];

	//	Vertex angles are in reversed order

	this.vA=new Array(3);
	for (var i=0;i<3;i++) this.vA[i]=neighbour.vA[2-i];

	//	Triangles are linked

	this.links=new Array(3);
	this.linkEdges=new Array(3);
	var ourEdge=2-edgeNumber;

	neighbour.links[edgeNumber]=this;
	neighbour.linkEdges[edgeNumber]=ourEdge;
	this.links[ourEdge]=neighbour;
	this.linkEdges[ourEdge]=edgeNumber;

	//	Two vertices are common

	var cv1=(edgeNumber+1)%3, cv2=(edgeNumber+2)%3;
	this.vertices=new Array(3);
	var c1=this.vertices[2-cv1]=neighbour.vertices[cv1];
	var c2=this.vertices[2-cv2]=neighbour.vertices[cv2];

	this.vertexOwner=new Array(3);
	this.vertexNumber=new Array(3);
	this.vertexOwner[ourEdge]=this;
	this.vertexNumber[ourEdge]=ourEdge;
	this.vertexOwner[2-cv1]=neighbour.vertexOwner[cv1];
	this.vertexNumber[2-cv1]=neighbour.vertexNumber[cv1];
	this.vertexOwner[2-cv2]=neighbour.vertexOwner[cv2];
	this.vertexNumber[2-cv2]=neighbour.vertexNumber[cv2];

	//	Third vertex is reflection in plane of shared edge.

	//	Unit normal to plane of reflection.

	var norm=neighbour.edgeNormals[edgeNumber];

	//	Construct the reflection in this unit normal; this is I - 2 n(x)n

	var reflection=Utils.multiDim([3,3],0);
	for (var i=0;i<3;i++)
	for (var j=0;j<3;j++)
		reflection[i][j]=(i==j?1.0:0.0)-2.0*norm[i]*norm[j];

	//	Reflect the neighbouring triangle's third vertex

	var c3=this.vertices[ourEdge]=new Array(3), v3=neighbour.vertices[edgeNumber];
	for (var i=0;i<3;i++)
		{
		var sum=0.0;
		for (var j=0;j<3;j++) sum+=reflection[i][j]*v3[j];
		c3[i]=sum;
		};
		
	//	Edge normals

	this.computeEdgeNormals();

	//	Form the isometry that takes this triangle back to the canonical position:
	//	this is the reflection, followed by the neighbour's isometry

	this.isometry=Utils.multiDim([3,3],0);
	var nIso=neighbour.isometry;
	for (var i=0;i<3;i++)
	for (var j=0;j<3;j++)
		{
		var sum=0.0;
		for (var k=0;k<3;k++) sum+=nIso[i][k]*reflection[k][j];
		this.isometry[i][j]=sum;
		};
	this.parity=-neighbour.parity;
	}
}

//	Create a tiling by reflection of this spherical triangle.

SphericalTriangle.prototype.tiling =
function (ntiles)
{
var tiles=new Array(ntiles);

tiles[0]=this;
this.tilingFlags=0;
this.nsteps=0;
this.steps=null;
var nt=1, firstFree=0;

while (nt<ntiles)
	{
	//	Find first free edge, and reflect in it to create another tile
	
	var tri, edge=0;
	findFree:
	for (tri=firstFree;tri<nt;tri++)
		for (edge=0;edge<3;edge++)
			if (tiles[tri].links[edge]==null) break findFree;
	if (tri==nt) break;
	
	firstFree=tri;
	var oldTile=tiles[tri];
	var newTile=tiles[nt]=new SphericalTriangle(oldTile,edge);
	oldTile.tilingFlags|=(1<<edge);
	
	var oldN=oldTile.nsteps, oldSteps=oldTile.steps;
	var newN=newTile.nsteps=oldN+1, newSteps=newTile.steps=new Array(newN);
	for (var i=0;i<oldN;i++) newSteps[i]=oldSteps[i];
	newSteps[oldN]=edge;
	
	var newVertex=2-edge;
	var nvtx=newTile.vertices, nvtx1=nvtx[newVertex];

	//	Check to see if the new vertex exists already
	
	scanVertices:
	for (var i=0;i<nt;i++)
		{
		oldTile=tiles[i];
		for (var j=0;j<3;j++)
			if (oldTile.vertexOwner[j]==oldTile
				&& SphericalTriangle.equals(nvtx1,oldTile.vertices[j]))
				{
				nvtx1=nvtx[newVertex]=oldTile.vertices[j];
				newTile.vertexOwner[newVertex]=oldTile;
				newTile.vertexNumber[newVertex]=j;
				break scanVertices;
				};
		};

	//	Check to see if either of the new tile's free edges match up with
	//	existing free edges.
	
	var correctParity=-newTile.parity;
	
	checkEdges:
	for (var nfree=1;nfree<=2;nfree++)
		{
		var newVertex2=(newVertex+nfree)%3;
		var nvtx2=nvtx[newVertex2];
		var newEdge=(newVertex+3-nfree)%3;
		var oldEdge=2-newEdge;
		
		for (tri=firstFree;tri<nt;tri++)
			{
			oldTile=tiles[tri];
			if (oldTile.parity==correctParity && oldTile.links[oldEdge]==null)
				{
				var vtx=oldTile.vertices;
				if (vtx[2-newVertex]==nvtx1 && vtx[2-newVertex2]==nvtx2)
					{
					newTile.links[newEdge]=oldTile;
					newTile.linkEdges[newEdge]=oldEdge;
					oldTile.links[oldEdge]=newTile;
					oldTile.linkEdges[oldEdge]=newEdge;
					continue checkEdges;
					};
				};
			};
		};
						
	nt++;
	};
return tiles;
};

//	Compute edge normals

SphericalTriangle.prototype.computeEdgeNormals =
function ()
{
this.edgeNormals=new Array(3);
for (var i=0;i<3;i++)
	{
	var j=(i+1)%3, k=(i+2)%3;
	this.edgeNormals[i]=Utils.cross(this.vertices[j],this.vertices[k]);
	Utils.normalise(this.edgeNormals[i]);
	};
}

//	Compare two unit vectors

SphericalTriangle.equals =
function (v, w)
{
var dp=-1.0;
for (var i=0;i<3;i++) dp+=v[i]*w[i];
return (Math.abs(dp)<1e-10);
}

//	Transform this triangle.

SphericalTriangle.prototype.rotate =
function (R)
{
for (var i=0;i<3;i++)
	{
	if (this.vertexOwner[i]==this) Utils.transform1(R, this.vertices[i]);
	Utils.transform1(R, this.edgeNormals[i]);
	};
	
//	Pre-multiply the isometry by R^{-1} = R^t
//	each original isometry S takes an original triangle to the canonical one,
//	so SR^{-1} takes a rotated triangle to the canonical one.

var rs=Utils.multiDim([3,3],0);
for (var i=0;i<3;i++)
for (var j=0;j<3;j++)
for (var k=0;k<3;k++)
	rs[i][j]+=this.isometry[i][k]*R[j][k];
this.isometry=rs;
}

//	inCentre()
//	==========

//	The point that lies on the bisectors of all the angles of the triangle.

SphericalTriangle.prototype.inCentre =
function ()
{
this.computeSideTrig();

//	Coefficients of incentre in terms of vertices

var denom=0.0, coeffs=new Array(3);
for (var i=0;i<3;i++)
	denom+=this.sideSin[i]*this.sideSin[i]+2.0*this.sideCos[i]*this.sideSin[(i+1)%3]*this.sideSin[(i+2)%3];
denom=Math.sqrt(denom);

for (var i=0;i<3;i++) coeffs[i]=this.sideSin[i]/denom;

//	Compute incentre

var c=[0,0,0];
for (var i=0;i<3;i++)
for (var j=0;j<3;j++)
	c[j]+=coeffs[i]*this.vertices[i][j];

return c;
}

//	pointFromAngles()
//	=================

//	Gives the coordinates of a point in the triangle from the three cosines of
//	the angles the point's vector makes with the vertices.

SphericalTriangle.prototype.pointFromAngles =
function (vcos)
{
this.computeSideTrig();

//	Point is p = x A + y B + z C, where A, B, C are vertices
//
//	Solve for x, y, z from values for p.A, p.B, p.C

var ca=this.sideCos[0], cb=this.sideCos[1], cc=this.sideCos[2];
var det=1.0-ca*ca-cb*cb-cc*cc+2.0*ca*cb*cc;
var mat=Utils.multiDim([3,3],0);
for (var i=0;i<3;i++)
for (var j=0;j<3;j++)
	{
	mat[i][j]=((i==j) ? this.sideSin[i]*this.sideSin[i] :
		this.sideCos[i]*this.sideCos[j]-this.sideCos[3-i-j]) / det;
	};
	
//	Coefficients of rotation point in terms of vertices

var coeffs=[0,0,0];
for (var i=0;i<3;i++)
for (var j=0;j<3;j++)
	coeffs[i]+=mat[i][j]*vcos[j];

//	Compute point

var pt=[0,0,0];
for (var i=0;i<3;i++)
for (var j=0;j<3;j++)
	pt[j]+=coeffs[i]*this.vertices[i][j];

return pt;
}

//	Compute cosines and sines of angular side lengths

SphericalTriangle.prototype.computeSideTrig =
function ()
{
this.sideCos=new Array(3);
this.sideSin=new Array(3);
for (var i=0;i<3;i++)
	{
	var dp=0, a=this.vertices[(i+1)%3], b=this.vertices[(i+2)%3];
	for (var l=0;l<3;l++) dp+=a[l]*b[l];
	this.sideCos[i]=dp;
	this.sideSin[i]=Math.sqrt(1.0-dp*dp);
	};
}
