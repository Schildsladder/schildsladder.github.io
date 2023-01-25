/*

Polyhedron.js
===============

Class representing a polyhedron.

*/

/* Original Java class members:

public int NF, NE, NV;				//	Number of faces, edges, vertices

public double vertices[][];			//	List of vertex (x,y,z) coordinates
public double faces[][];			//	List of face outward normals
									//	4th component in each array is distance of face
									//	from origin.

public int face_vertices[][];		//	List of face vertex numbers
public int face_edges[][];			//	List of face edge numbers
public boolean face_edge_rev[][];	//	True if edge points backwards compared to
									//	ascending face vertex numbers

public int edge_vertices[][];		//	Vertex numbers for each edge
public int edge_faces[][];			//	Face numbers for each edge

public int vertex_vertices[][];		//	List of vertices around each vertex
public int vertex_edges[][];		//	List of edges around each vertex
public int vertex_faces[][];		//	List of faces around each vertex

public double isometries[][][];
public boolean even_isometry[];
public int vertex_permutations[][];
public int face_permutations[][];
*/


//	Constructors
//	============

//	(1.) Construct a polyhedron from an array of vertices.
//
//	vertices			gives the vertex coordinates
//
//	nv					is the count of vertices
//
//	vertex_vertices		indicates which vertex is joined by an edge to which.
//						The list for each vertex must be ordered counterclockwise,
//						and consecutive entries must have a face in common

function Polyhedron(vertices, nv, vertex_vertices)
{
this.vertices=vertices;					//	Store the vertex coordinates,
this.NV=nv;									//	vertex count,
this.vertex_vertices=vertex_vertices;	//	and list of vertex-vertex connections.

this.vertex_edges=new Array(nv);				//	Arrays with details for each vertex
this.vertex_faces=new Array(nv);
var ne=0;								//	Count total number of edges
for (var i=0;i<nv;i++)
	{
	var nvv=vertex_vertices[i].length;
	this.vertex_edges[i]=new Array(nvv);
	this.vertex_faces[i]=new Array(nvv);
	ne+=nvv;
	};
this.NE=ne=ne/2;								//	Every edge was counted twice

this.edge_vertices=Utils.multiDim([ne,2],0);			//	Arrays with details for each edge
this.edge_faces=Utils.multiDim([ne,2],0);

//	Generate edges

var evv=Utils.multiDim([nv,nv],0);			//	Edge for each vertex pair
var vvn=Utils.multiDim([nv,nv],0);			//	Link number for each vertex wrt each
var je=0;								//	Count edges as they're created
for (var i=0;i<nv;i++)
	{
	var vvi=vertex_vertices[i], nvv=vvi.length;
	for (var j=0;j<nvv;j++)
		{
		var k=vvi[j];					//	Vertex i and k are linked
		vvn[i][k]=j;					//	Log position at which vertex k connects to i
		if (k>i)						//	Second vertex is greater, create edge
			{
			this.edge_vertices[je][0]=i;
			this.edge_vertices[je][1]=k;
			evv[i][k]=evv[k][i]=je;
			je++;
			};
		};
	};

//	Record the edges around each vertex

for (var i=0;i<nv;i++)
	{
	var vvi=vertex_vertices[i], nvv=vvi.length;
	for (var j=0;j<nvv;j++) this.vertex_edges[i][j]=evv[i][vvi[j]];
	};

//	Generate faces

var nf=ne;								//	Maximum possible number of faces
this.face_vertices=new Array(nf);			//	Arrays with details for each face
this.face_edges=new Array(nf);
this.face_edge_rev=new Array(nf);
this.faces=Utils.multiDim([nf,4],0);

//	Mark vertex face slots as all empty

for (var i=0;i<nv;i++)
	{
	var mf=this.vertex_faces[i].length;
	for (var j=0;j<mf;j++) this.vertex_faces[i][j]=-1;
	};
	
var jf=0;
var fv=new Array(nv), fe=new Array(nv);
var rev=new Array(nv);
for (var tt1=0;tt1<10000;tt1++)
	{
	//	Find vertex, i, with unused starting edge, k
	
	var i=0, k=0;
	vLoop:
	for (i=0;i<nv;i++)
		{
		var mf=this.vertex_faces[i].length;
		for (k=0;k<mf;k++)
			if (this.vertex_faces[i][k]<0) break vLoop;
		};
	if (i==nv) break;
	
	var ns=0;
	for (var tt2=0;tt2<10000;tt2++)
		{
		fv[ns]=i;							//	Log vertex as part of face
		this.vertex_faces[i][k]=jf;				//	Log face as adjoining vertex
		var m=vertex_vertices[i][k];		//	Find next vertex
		var n=evv[i][m];
		fe[ns]=n;							//	Log edge
		rev[ns]=(i != this.edge_vertices[n][0]);	//	First vertex not first in edge
		this.edge_faces[n][rev[ns]?1:0]=jf;		//	Log face as adjoining this edge
		ns++;
		
		var mf=vertex_vertices[m].length;
		k=(vvn[m][i]+mf-1)%mf;				//	Move clockwise around new vertex
		i=m;								//	And go to it
		if (m==fv[0]) break;				//	Unless we've already been there
		};

	//	Compute face normals/distances from origin

	var fn=this.faces[jf];
	
	//	Three consecutive corners of face are a, b, c
	
	var a=this.vertices[fv[0]], b=this.vertices[fv[1]], c=this.vertices[fv[2]];
	
	//	Compute normal to face as (b-a)X(c-b)
	
	for (var j=0;j<3;j++)
		{
		var l=(j+1)%3, m=(j+2)%3;
		fn[j]=(b[l]-a[l])*(c[m]-b[m])-(b[m]-a[m])*(c[l]-b[l]);
		};
	Utils.normalise(fn);
	
	//	Dot product of normal with one corner is distance from origin
	
	fn[3]=Utils.dot(fn,a);
	
	var fvjf=this.face_vertices[jf]=new Array(ns);
	var fejf=this.face_edges[jf]=new Array(ns);
	var revjf=this.face_edge_rev[jf]=new Array(ns);
	for (var j=0;j<ns;j++)
		{
		fvjf[j]=fv[j];
		fejf[j]=fe[j];
		revjf[j]=rev[j];
		};
	jf++;
	};
this.NF=nf=jf;

this.isometries=null;
this.even_isometry=null;
this.vertex_permutations=null;
this.face_permutations=null;
}

//	

//	Create a convex regular polyhedron with vertices on the unit sphere.
//	One vertex is at (0,0,1) and one edge from that vertex is aligned in the +ve
//	x-direction (as well as being tilted to some degree in the z-direction).

//	Wythoff symbol data for the five regular polyhedra

Polyhedron.regularWythoff=[
[[1,3],[1,2],[1,3],[1,24]],		//	tetrahedron
[[1,3],[1,2],[1,4],[1,48]],		//	cube
[[1,4],[1,2],[1,3],[1,48]],		//	octahedron
[[1,3],[1,2],[1,5],[1,120]],	//	dodecahedron
[[1,5],[1,2],[1,3],[1,120]]		//	icosahedron
];

Polyhedron.regularPolyhedron =
function (index)
{
return Polyhedron.wythoffPolyhedron(Polyhedron.regularWythoff[index], false, false, false);
};

//	Create a Polyhedron from its Wythoff symbol

Polyhedron.wythoffPolyhedron =
function (params,
	storeIsometries,
	vertexPermutations,
	facePermutations)
{
var vA=params, ntiles=params[3][1], wythoffKind=params[3][0];
var n=vA[0][1], d=vA[0][0];

var canonical=new SphericalTriangle(vA);
var tiles=canonical.tiling(ntiles);

var vertices=null;
var vv=null, nv=0, mf;

//	| p q r

if (wythoffKind==0)
	{
	//	Vertices of polyhedron lie at the point that, when rotated by twice the
	//	vertex angle around any of the vertices, is displaced by the same angle in
	//	each case -- which is the angle subtended by an edge of the polyhedron.
	
	//	Work out vertex configuration
	
	mf=3;
	var rightAngle=[false,false,false];
	for (var k=0;k<3;k++)
		if (!(rightAngle[k]=(vA[k][1]==2))) mf++;
		
	var ns=new Array(mf);
	var nsi=0;
	for (var k=0;k<3;k++)
		{
		ns[nsi++]=3.0;
		if (!rightAngle[k]) ns[nsi++]=vA[k][1]/vA[k][0];
		};
		
	//	Find the half edge angle, then cos of full edge angle
	
	var hea=Polyhedron.computeGeometry(ns,1,null,null);
	var cosea=Math.cos(2.0*hea);
	
	//	Find the cosines of the angles separating the polyhedron vertex from
	//	the three triangle vertices, given the known displacement angle when
	//	rotated by twice the angles at the vertices.
	
	var ca=new Array(3);
	for (var k=0;k<3;k++)
		{
		var cRot=Math.cos(2.0*Math.PI*vA[k][0]/vA[k][1]);
		ca[k]=Math.sqrt((cosea-cRot)/(1.0-cRot));
		};
		
	//	Find the point from these angles

	var canonicalVertex=canonical.pointFromAngles(ca);
	
	//	Use isometries to generate the other vertices, in all triangles with
	//	parity of 1.

	nv=Math.floor(ntiles/2);
	vertices=Utils.multiDim([nv,3],0);
	var vtiles=new Array(nv);
	var jv=0;
	for (var i=0;i<ntiles;i++)
		{
		var t=tiles[i];
		if (t.id==0 && t.parity==1)
			{
			vtiles[jv]=t;
			for (var j=0;j<3;j++)
				{
				var sum=0.0;
				for (var k=0;k<3;k++) sum+=t.isometry[k][j]*canonicalVertex[k];
				vertices[jv][j]=sum;
				};
			t.id=1+jv;
			jv++;
			};
		};
		
	//	Identify connections between vertices, twice reflected

	vv=Utils.multiDim([nv,mf],0);
	for (var i=0;i<nv;i++)
		{
		var t=vtiles[i];	//	Start with triangle that owns vertex
		var j=0;
		for (var k=0;k<3;k++)
			{
			var u=t.links[k];
			vv[i][j++]=u.links[(2-k+1)%3].id-1;
			if (!rightAngle[(k+2)%3]) vv[i][j++]=u.links[(2-k+2)%3].id-1;
			};
		};
	}
	
//	p | q r

else if (wythoffKind==1)
	{
	//	Vertices of polyhedron lie at point P of triangle PQR
	
	nv=Math.floor(ntiles/(2*n));
	vertices=new Array(nv);
	var vtiles=new Array(nv);
	var jv=0;
	for (var i=0;i<ntiles;i++)
		{
		var t=tiles[i];
		if (t.id==0)
			{
			vtiles[jv]=t;
			vertices[jv]=t.vertices[t.parity==1?0:2];

			//	Give all triangles that share vertex P the same ID
			
			var u=t;
			do
				{
				u.id=1+jv;
				}
			while ((u=u.links[u.parity==1?1:0])!=t);
			jv++;
			};
		};
		
	//	Identify connections between vertices, reflected in side opposite P

	var qRightAngle=(vA[1][1]==2);
	mf = qRightAngle ? n : 2*n;
	vv=Utils.multiDim([nv,mf],0);
	for (var i=0;i<nv;i++)
		{
		var t=vtiles[i];	//	Start with triangle that owns vertex
		for (var j=0;j<mf;j++)			//	Move around vertex, counterclockwise
			{
			var k,l,m;
			if (t.parity==1) {k=0; l=1; m=0;} else {k=2; l=0; m=1;};
			vv[i][j]=t.links[k].id-1;
			t=t.links[l];
			if (qRightAngle) t=t.links[m];
			};
		};
	}

//	p q | r

else if (wythoffKind==2)
	{
	//	Vertices of polyhedron lie at point on PQ that bisects angle at R
	
	//	Vertex angles, halving the third one

	var a=Math.PI*vA[0][0]/vA[0][1];
	var c=Math.PI*vA[2][0]/vA[2][1]/2.0;

	//	Sines and cosines

	var ca=Math.cos(a), sa=Math.sin(a);
	var cc=Math.cos(c), sc=Math.sin(c);
	
	//	Geodesic angles
	
	var cosBeta=canonical.vertices[2][2];
	var sinBeta=Math.sqrt(1.0-cosBeta*cosBeta);
	
	var gamma=Math.atan(sinBeta/(ca*cosBeta+sa*cc/sc));
	
	//	Canonical vertex

	var canonicalVertex=[Math.sin(gamma),0.0,Math.cos(gamma)];
	
	//	Use isometries to generate the other vertices
	
	nv=Math.floor(ntiles/2);
	vertices=Utils.multiDim([nv,3],0);
	var vtiles=new Array(nv);
	var jv=0;
	for (var i=0;i<ntiles;i++)
		{
		var t=tiles[i];
		if (t.id==0 && t.parity==1)
			{
			vtiles[jv]=t;
			for (var j=0;j<3;j++)
				{
				var sum=0.0;
				for (var k=0;k<3;k++) sum+=t.isometry[k][j]*canonicalVertex[k];
				vertices[jv][j]=sum;
				};
				
			//	Give the other triangle, reflection in PQ, the same ID
			
			t.links[2].id=t.id=1+jv;
			jv++;
			};
		};
		
	//	Identify connections between vertices, reflected in sides other than PQ

	var pRightAngle=(vA[0][1]==2);
	mf = pRightAngle ? 3 : 4;
	vv=Utils.multiDim([nv,mf],0);
	for (var i=0;i<nv;i++)
		{
		var t=vtiles[i];	//	Start with triangle that owns vertex
		vv[i][0]=t.links[0].id-1;
		vv[i][1]=t.links[1].id-1;
		if (!pRightAngle) vv[i][2]=t.links[2].links[1].id-1;
		vv[i][mf-1]=t.links[2].links[2].id-1;
		};
	}
	
//	p q r |

else if (wythoffKind==3)
	{
	//	Vertices of polyhedron lie at incentre of PQR
	
	var canonicalVertex=canonical.inCentre();
	
	nv=ntiles;
	vertices=Utils.multiDim([nv,3],0);
	for (var i=0;i<ntiles;i++)
		{
		var t=tiles[i];
		for (var j=0;j<3;j++)
			{
			var sum=0.0;
			for (var k=0;k<3;k++) sum+=t.isometry[k][j]*canonicalVertex[k];
			vertices[i][j]=sum;
			};
		t.id=1+i;
		};
		
	//	Identify connections between vertices

	mf=3;
	vv=Utils.multiDim([nv,mf],0);
	for (var i=0;i<nv;i++)
		{
		var t=tiles[i];
		for (var k=0;k<3;k++)
			vv[i][k]=t.links[k].id-1;
		};
	};

var ph=new Polyhedron(vertices, nv, vv);

//	Retain the isometries in the Polyhedron object.

var chiral=(wythoffKind==0);				//	Only uses even parity tiles
var nIsometries=chiral ? Math.floor(ntiles/2) : ntiles;

if (storeIsometries)
	{
	ph.isometries=new Array(nIsometries);
	ph.even_isometry=new Array(nIsometries);
	var nIso=0;
	for (var i=0;i<ntiles;i++)
		{
		var even=(tiles[i].parity==1);
		if (!chiral || even)
			{
			ph.isometries[nIso]=tiles[i].isometry;
			ph.even_isometry[nIso++]=even;
			};
		};
	};
	
//	Generate a list of permutations of the vertices, corresponding
//	to each isometry.

if (vertexPermutations || facePermutations)
	{
	ph.vertex_permutations=Utils.multiDim([nIsometries,nv],0);
	
	//	Visit all tiles from the canonical tile, and shadow all these
	//	moves from the isometry tile.  The result maps vertices to vertices.
	
	var nIso=0;
	for (var i=0;i<ntiles;i++)
		if (!chiral || tiles[i].parity==1)
			Polyhedron.visitTile(canonical, tiles[i], tiles[i].parity, 
				ph.vertex_permutations[nIso++], chiral);
	};
	
//	Generate a list of permutations of the faces, corresponding
//	to each isometry.

if (facePermutations)
	{
	var nf=ph.NF;
	ph.face_permutations=Utils.multiDim([nIsometries,nf],0);
	
	//	Create a bit-per-vertex signature for each face

	var fsn=1+Math.floor(nv/31);
	var face_signatures=new Array(nf);
		
	//	For each isometry, use the vertex permutations to deduce face permutations
	
	for (var nIso=0;nIso<nIsometries;nIso++)
		{
		var vperm=ph.vertex_permutations[nIso];
		var fperm=ph.face_permutations[nIso];
		for (var i=0;i<nf;i++)
			{
			var fvi=ph.face_vertices[i],nfv=fvi.length, psig=new Array(fsn);
			for (var j=0;j<nfv;j++)
				{
				var vn=vperm[fvi[j]];
				psig[Math.floor(vn/31)]|=(1<<(vn%31));
				};
			if (nIso==0)
				{
				face_signatures[i]=psig;
				fperm[i]=i;
				}
			else for (var fp=0;fp<nf;fp++)
				{
				var k;
				for (k=0;k<fsn;k++) if (psig[k]!=face_signatures[fp][k]) break;
				if (k==fsn)
					{
					fperm[i]=fp;
					break;
					};
				};
			};
		};
	if (!vertexPermutations) ph.vertex_permutations=null;
	};

return ph;
}

Polyhedron.visitTile =
function (t, shadow, parity, perm, chiral)
{
if (!chiral || t.parity==1) perm[t.id-1]=shadow.id-1;
for (var ie=0;ie<3;ie++)
	if ((t.tilingFlags & (1<<ie))!=0)
		Polyhedron.visitTile(t.links[ie], shadow.links[parity==1?ie:2-ie], parity, perm, chiral);
}

//	computeGeometry()
//	-----------------

//	Compute the geometry of a uniform polyhedron with a specified vertex configuration
//	and vertex density.
//	
//	This method is based on "Uniform Solution for Uniform Polyhedra" by Zvi Har'El.
//
//	Terminology
//	===========
//
//	Each face is a regular polygon, with some quantities defined as follows:
//
//	- ns, the number of sides.  This can be modified for some special cases:
//		* For star faces, ns->ns/(face density), where the face density is the
//		number of 2pi rotations around the face's centre that the ns sides comprise.
//		* For retrograde faces, ns->ns/(ns-1).  This converts alpha=PI/ns,
//		to PI-alpha, and reverses the sign of gamma.
//
//	- alpha, half the angle that each side of the face subtends at the face's centre.
//	- gamma, half the angle that each corner of the face occupies at the vertex.
//	- A, half the angular geodesic length of each side (common to all faces)
//	- B, the angular geodesic distance from the vertex to the centre of the face.
//
//	The vertex configuration is supplied as an array, ns, giving the number of sides
//	in the sequence of regular polygonal faces found around every vertex (possibly
//	modified as described above).
//
//	The vertex density supplied is an integer, vDensity, which gives the number of
//	rotations around the vertex that the sequence of polygons comprises.
//
//	Returns:
//	--------
//
//	The function value returned is A, half the angular edge length.
//
//	The array faceGamma is set, for each face in the supplied sequence, to half
//	the angle that a corner of the face occupies at the vertex.
//
//	The array faceB is set, for each face in the supplied sequence, to the angular
//	length of the geodesic from the vertex to the centre of the face.

Polyhedron.computeGeometry =
function (ns, vDensity, faceGamma, faceB)
{
var tolerance=1e-10;				//	Required accuracy

var m=ns.length;					//	Number of faces around the vertex
var nsD=Utils.multiDim([m],0);				//	Distinct values from ns
var mult=Utils.multiDim([m],0);				//	Multiplicity for each entry in nsD
var numD=0;							//	Count of these distinct values
var maxNS=0;						//	Index to maximum value in nsD
var faceKind=new Array(m);			//	Index to nsD for each face in sequence

for (var i=0;i<m;i++)				//	For each face in sequence
	{
	var nsi=ns[i];				//	Get its number of sides
	var j;
	for (j=0;j<numD;j++)
		if (nsD[j]==nsi) break;		//	See if it's already in the list, nsD
	if (j==numD)
		{
		nsD[numD++]=nsi;			//	It's not, so add it
		if (nsi>nsD[maxNS])
			maxNS=j;				//	Track maximum
		};
	faceKind[i]=j;					//	Either way, record the index
	mult[j]++;						//	and update its multiplicity
	};
	
var cosAlpha=new Array(numD);	//	Alpha is half the angle a face's side subtends
									//	at its centre
									
var gamma=new Array(numD);		//	Gamma is half the angle a face's corner occupies
									//	at the vertex
									
var gammaSum=Math.PI*vDensity;	//	Multiplicity-weighted sum of gammas
									
var cosA;						//	A is half the angular length of every edge

for (var i=0;i<numD;i++)
	{
	var alpha=Math.PI/nsD[i];
	cosAlpha[i]=Math.cos(alpha);
	};

//	If there is only one kind of face, we can immediately solve.

if (numD==1)
	{
	gamma[0]=gammaSum/m;
	cosA=cosAlpha[0]/Math.sin(gamma[0]);
	}
	
//	There is more than one kind of face, so we find a solution iteratively.

else
	{
	//	Set up initial values for gammas, based on Euclidean triangles with angles
	//	summing to pi.
	
	for (var i=0;i<numD;i++) gamma[i]=Math.PI*(0.5-1.0/nsD[i]);
	cosA=cosAlpha[maxNS]/Math.sin(gamma[maxNS]);
	
	//	Solve iteratively
	
	for (var tt3=0;tt3<10000;tt3++)
		{
	
		//	Compute delta, the difference between target sum and actual sum of gammas
	
		var delta=gammaSum;
		for (var i=0;i<numD;i++) delta-=mult[i]*gamma[i];
		
		//	If the difference is small enough, we're finished
		
		if (Math.abs(delta)<tolerance) break;
		
		//	Update the independent variable gamma[maxNS]
		
		var sumTan=0.0;
		for (var i=0;i<numD;i++) sumTan+=mult[i]*Math.tan(gamma[i]);
		gamma[maxNS]+=delta*Math.tan(gamma[maxNS])/sumTan;
		
		//	Update the dependent variables
		
		cosA=cosAlpha[maxNS]/Math.sin(gamma[maxNS]);
		for (var i=0;i<numD;i++)
			if (i!=maxNS) gamma[i]=Math.asin(cosAlpha[i]/cosA);
		};
	};
	
//	Compute the angles B from the vertex to the centre of the faces

var B=new Array(numD);
for (var i=0;i<numD;i++)
	B[i]=Math.acos(1.0/(
		Math.tan(Math.PI/nsD[i])*Math.tan(gamma[i])
			));
	
//	Set return values for all faces
	
for (var i=0;i<m;i++)
	{
	var k=faceKind[i];
	if (faceGamma!=null) faceGamma[i]=gamma[k];
	if (faceB!=null) faceB[i]=B[k];
	};

return Math.acos(cosA);
}

//	Transform the polyhedron with a given rotation matrix and scale factor
//	======================================================================

Polyhedron.prototype.transform =
function (R, scale)
{
var tmp=new Array(3);
var tl=[this.vertices, this.faces];	//	Arrays we need to transform
var tn=[this.NV,this.NF];
var reScale=(scale!=1.0);

for (var vf=0;vf<2;vf++)
	{
	var reScaleThis=reScale && vf==0;
	var array=tl[vf];
	var na=tn[vf];
	for (var i=0;i<na;i++)
		{
		var ai=array[i];
		for (var j=0;j<3;j++)
			{
			tmp[j]=0.0;
			var Rj=R[j];
			for (var k=0;k<3;k++) tmp[j]+=Rj[k]*ai[k];
			};
		for (var j=0;j<3;j++) ai[j]=(reScaleThis ? scale*tmp[j] : tmp[j]);
		};
	};
if (reScale) for (var i=0;i<this.NF;i++) this.faces[i][3]*=scale;
}

//	Decide if a given point is in the interior of this polyhedron
//	=============================================================

Polyhedron.prototype.contains =
function (pt)
{
var nf=this.NF;
for (var i=0;i<nf;i++)		//	For all faces
	{
	//	Form dot product of face normal with point vector,
	//	minus distance of face from origin.
	
	var fn=this.faces[i];
	var dp=-fn[3];
	for (var k=0;k<3;k++) dp+=pt[k]*fn[k];
	if (dp>0.0) return false;
	};
return true;
}

//	Find all the intersections between this polyhedron's edges and another polyhedron
//	=================================================================================

//	Returns the number of intersections.
//
//	Count of intersections for edge ie, followed by face number/point number pairs,
//	are placed in eIdata[ie][0-5]
//	Pointers to the intersection points are placed in pts

Polyhedron.prototype.edgeIntersections =
function (phB, eIdata, pts)
{
//	Compute distances from all phB faces to all of this polyhedron's vertices

var nv=this.NV, nfB=phB.NF;
var fvd=Utils.multiDim([nv,nfB],0);
for (var iv=0;iv<nv;iv++)
	{
	var vv=this.vertices[iv], fv=fvd[iv];
	for (var f=0;f<nfB;f++)
		{
		var fn=phB.faces[f], dp=-fn[3];
		for (var k=0;k<3;k++) dp+=fn[k]*vv[k];
		fv[f]=dp;
		};
	};

//	Find intersections along edges

var ee=new Array(2);
var ne=this.NE, fi=new Array(2), nInter=0;

edgeLoop:
for (var ie=0;ie<ne;ie++)
	{
	var eIe=eIdata[ie];
	eIe[0]=0;
	var ev=this.edge_vertices[ie], vi0=ev[0], vi1=ev[1];
	var fvd0=fvd[vi0], fvd1=fvd[vi1];
	var entry=0.0, exit=1.0;
	var fin=0, fout=0;
	
	for (var f=0;f<nfB;f++)		//	For all faces of phB
		{
		var dpp=fvd0[f], dpd=fvd1[f]-dpp;
			
		if (dpd==0.0)	//	Edge parallel to this face
			{
			if (dpp>0.0)
				continue edgeLoop;		//	Starting point is on outside of this face-plane
			}
		else
			{
			var lambda=-dpp/dpd;		//	Value at which v0 + lambda (v1-v0) hits face-plane
			
			if (dpd<0.0)				//	Edge is oriented going INTO face
				{
				if (lambda>entry)		//	Later entry than we have from other faces?
					{entry=lambda; fin=f;};
				}
			else						//	Edge is oriented going OUT of face
				{
				if (lambda<exit)		//	Earlier exit than we have from other faces?
					{exit=lambda; fout=f;};
				};
			if (entry>exit)
				continue edgeLoop;		//	No overlap of lambda left inside polyhedron
			};
		};
		
	ee[0]=entry;
	ee[1]=exit;
	fi[0]=fin;
	fi[1]=fout;
	var v0=this.vertices[vi0], v1=this.vertices[vi1];
	var nei=0;
	for (var j=0;j<2;j++)
		{
		var lambda=ee[j];
		if (lambda>0.0 && lambda<1.0)
			{
			var newPt=new Array(3), oml=1.0-lambda;
			for (var k=0;k<3;k++) newPt[k]=oml*v0[k]+lambda*v1[k];
			pts[nInter]=newPt;
			eIe[2*nei+1]=fi[j];
			eIe[2*(nei++)+2]=nInter++;
			};
		};
	eIe[0]=nei;
	};
return nInter;
}

//	Split the faces of two polyhedra into parts lying inside/outside each other
//	===========================================================================
//
//	Returns a Vector containing Gon3D objects, which describe the polygons created
//	from splitting the faces.
//
//	The "surfaceID" of each Gon3D has the face number of its parent face.
//
//	The "flags" of each Gon3D are set from:
//
//		POLYHEDRON_0			a piece of a face of the first polyhedron
//		POLYHEDRON_1			a piece of a face of the second polyhedron
//		INTERIOR_PIECE			a piece that lies inside the other polyhedron
//		EXTERIOR_PIECE			a piece that lies outside the other polyhedron
//		FRONT_FACE				a piece of a face whose normal has a +ve z-component
//		BACK_FACE				a piece of a face whose normal has a 0/-ve z-component
//
//	The "edgeID" has flags set:
//
//		GonND.EDGE_ORIGINAL		all or part of one of the polyhedron's original edges
//		GonND.EDGE_INTERNAL		an edge created to split an exterior piece with a hole
//								in the middle into two polygons
//		GonND.EDGE_POSITIVE		OR'd with number of face that created edge
//		GonND.EDGE_NEGATIVE		ditto
//

Polyhedron.POLYHEDRON_0 = 		1<<0;
Polyhedron.POLYHEDRON_1 = 		1<<1;
Polyhedron.ALL_POLYHEDRA =		Polyhedron.POLYHEDRON_0 | Polyhedron.POLYHEDRON_1;
Polyhedron.INTERIOR_PIECE = 	1<<2;
Polyhedron.EXTERIOR_PIECE = 	1<<3;
Polyhedron.ALL_PIECES =			Polyhedron.INTERIOR_PIECE | Polyhedron.EXTERIOR_PIECE;
Polyhedron.FRONT_FACE = 		1<<4;
Polyhedron.BACK_FACE = 			1<<5;
Polyhedron.ALL_FACES =			Polyhedron.FRONT_FACE | Polyhedron.BACK_FACE;

Polyhedron.faceIntersections =
function (
	ph,							//	Array of two polyhedra
	choiceFlags					//	Flags saying which pieces we want
	)
{
var wantInterior=((choiceFlags & Polyhedron.INTERIOR_PIECE)!=0);
var wantExterior=((choiceFlags & Polyhedron.EXTERIOR_PIECE)!=0);

//	Original edgeIntersection() data:

var eIdataP=new Array(2);		//	Intersection data for the two polyhedra
var ptsP=new Array(2);	//	Intersection points for the two polyhedra

//	Get edgeIntersection() data, and put it into an array indexed by first/second faces.

//	First number is count of data, then there are up to two successive entries of:
//		- polyhedron number    }
//		- edge number          } for edge on which intersection point lies
//		- intersection number  }
//		- point number         }
//		- other face number associated with this edge
//		- index to other data entry associated with this edge

var PH=1, EDGE=2, INTER=3, PT=4, FACE=5, INDX=6, NDATA=6;
var nf1=ph[0].NF, nf2=ph[1].NF;

//	Index this by either face first, for convenience.

var idataP=new Array(2);
var idata0=idataP[0]=Utils.multiDim([nf1,nf2,2*NDATA+1],0);
var idata1=idataP[1]=Utils.multiDim([nf2,nf1],null);
for (var f1=0;f1<nf1;f1++) for (var f2=0;f2<nf2;f2++) idata1[f2][f1]=idata0[f1][f2];

//	Count edge intersections for each face

var eif=new Array(2);
	
for (var ip=0;ip<2;ip++)
	{
	var p=ph[ip];
	var ne=p.NE, eIdata=eIdataP[ip]=Utils.multiDim([ne,5],0), idata=idataP[ip];
	ptsP[ip]=[];
	p.edgeIntersections(ph[1-ip], eIdata, ptsP[ip]);
	var eifi=eif[ip]=Utils.multiDim([p.NF],0);
	
	for (var ie=0;ie<ne;ie++)			//	For each edge
		{
		var eId=eIdata[ie], ni=eId[0];
		var ef=p.edge_faces[ie];		//	List of faces for this edge
		var f1a=ef[0], f1b=ef[1];		//	The two faces of p that share the edge

		for (var n=0;n<ni;n++)
			{
			eifi[f1a]++;						//	Count edge intersections for each face
			eifi[f1b]++;
			var f2=eId[2*n+1];					//	Face number of other polyhedron
			var ida=idata[f1a][f2], idb=idata[f1b][f2];
			var da=ida[0], db=idb[0];			//	Indices to data
			ida[da+PH]=idb[db+PH]=ip;			//	Log polyhedron owning edge
			ida[da+EDGE]=idb[db+EDGE]=ie;		//	Log edge number
			ida[da+INTER]=idb[db+INTER]=n;		//	Log intersection number
			ida[da+PT]=idb[db+PT]=eId[2*n+2];	//	Log point number
			ida[da+FACE]=f1b;	idb[db+FACE]=f1a;	//	Log other face number
			ida[da+INDX]=db;	idb[db+INDX]=da;	//	Log other data index
			ida[0]+=NDATA; idb[0]+=NDATA;
			};
		};
	};

var gons=[];
	
for (var ip=0;ip<2;ip++)			//	For both polyhedra
	{
	var phFlag=1<<ip;				//	Flag that codes for this polyhedron
	if ((phFlag & choiceFlags)==0) continue;
	
	var p=ph[ip], q=ph[1-ip];
	var nf=p.NF, nfq=q.NF, maxgv=p.NV+2*nfq;
	var vertices=p.vertices;
	var pts=ptsP[ip];
	var gf=new Array(maxgv);
	var eids=new Array(maxgv);
	var eIdata=eIdataP[ip], idata=idataP[ip];
	
	for (var f1=0;f1<nf;f1++)			//	For all faces f1 of current polyhedron
		{
		var ns=p.face_vertices[f1].length;
		
		//	Skip if face is unwanted for being front/back
		
		var faceFlag=p.faces[f1][2]>0 ? Polyhedron.FRONT_FACE : Polyhedron.BACK_FACE;
		if ((faceFlag & choiceFlags)==0) continue;

		var idF=idata[f1];			//	Row of intersection data for this face
		var fv=p.face_vertices[f1], fe=p.face_edges[f1];
		
		//	Find first face of other polyhedron, ff, that cuts this face, f1
		
		var ff;
		for (ff=0;ff<nfq;ff++) if (idF[ff][0]==2*NDATA) break;
		
		if (ff<nfq)
			{
			
			var noEdgeIntersections=(eif[ip][f1]==0);
			var interiorPiece=null;
			
			//	Flag all polyhedron vertices for this face as unused
			
			var used=Utils.multiDim([ns],false);
	
			for (var pass=0;pass<2;pass++)	//	Interior/exterior
				{
				//	Skip interior or exterior pieces if they're not wanted;
				//	but if there are no edge intersections, we force calculation
				//	of interior piece, and use that to split
				//	exterior.

				var pieceFlag = pass==0 ? Polyhedron.INTERIOR_PIECE : Polyhedron.EXTERIOR_PIECE;
				if (!(noEdgeIntersections && pass==0) &&
					(pieceFlag & choiceFlags)==0) continue;
				var edgeFlag = pass==0 ? GonND.EDGE_NEGATIVE : GonND.EDGE_POSITIVE;
			
				//	Exterior of face with no edge intersections has to be returned
				//	as two artificial pieces, because it's a polygon with a hole
				//	in the middle.
					
				if (pass==1 && noEdgeIntersections)
					{
					//	Details of interior piece

					var nip=interiorPiece.nvert, ipID=interiorPiece.edgeIDs;
					var ipV=interiorPiece.vertices;
					
					//	Split this face on first face of other polyhedron that intersects

					var gnd=p.faceToGon3D(f1);
					var vNeg=[], vOn=[], vPos=[];
					gnd.BSPsplit(q.faces[ff],GonND.EDGE_INTERNAL,null,0,null,vNeg,vOn,vPos);
					
					//	Append appropriate parts of perimeter of interior piece to
					//	the two halves
					
					for (var aa=0;aa<2;aa++)
						{
						var vvaa=(aa==0?vNeg:vPos);
						if (vvaa.length==1)
							{
							var half=vvaa[0];
							var hV=half.vertices;
							var nh=half.nvert, hID=half.edgeIDs, nap, v1, v2;
							if (aa==1) {nap=nh+2; v1=1; v2=0;}
							else {nap=nh+nip; v1=0; v2=1;}
							
							var napV=new Array(nap);
							var napID=new Array(nap);
							
							var napC=0;
							for (var i=0;i<nh;i++)
								{
								napV[napC]=hV[i];
								var hi=hID[i];
								napID[napC++]=hi;
								if ((hi & GonND.EDGE_INTERNAL)!=0)
									{
									var vi=v1;
									for (var tt3=0;tt3<10000;tt3++)
										{
										napV[napC]=ipV[vi];
										napID[napC++]=ipID[(vi+nip-1)%nip];
										if (vi==v2) break;
										vi=(vi+nip-1)%nip;
										};
									napID[napC-1]=GonND.EDGE_INTERNAL;
									};
								};
							var mpiece=new Gon3D(napV,nap);
							mpiece.edgeIDs=napID;
							mpiece.surfaceID=f1;
							mpiece.flags=phFlag | pieceFlag | faceFlag;
							gons.push(mpiece);
							};
						};
					}
				
				//	Ordinary cases
				
				else for (var tt4=0;tt4<10000;tt4++)						//	For individual pieces
					{
					var indx=0, ngf=0, vi=0, id=null, ie=0, ien=0, f2=0;
					var p0=null, pn=null;
					var needPoint, onInter;
					
					if (pass==0)
						{
						//	Interior piece; start with segment where face ff cuts face f1
						
						//	Use scalar triple product with face normals to get direction
						//	to traverse the line segment where ff cuts f1, such that
						//	the interior of the other polyhedron lies on the left.
						
						id=idF[ff];
						var fc1=p.faces[f1], fc2=q.faces[ff];
						var pt1=ptsP[id[PH]][id[PT]];
						var pt2=ptsP[id[NDATA+PH]][id[NDATA+PT]];
						var stp=0.0;
						for (var i=0;i<3;i++)
							{
							var j=(i+1)%3;
							var k=(i+2)%3;
							stp+=(pt2[i]-pt1[i])*(fc1[j]*fc2[k]-fc1[k]*fc2[j]);
							};

						if (stp>0) {indx=0; p0=pt1;} else {indx=NDATA; p0=pt2;};
						needPoint=false;
						onInter=true;
						f2=ff;
						}
					else
						{
						//	Exterior piece; start with first unused vertex;
						//	if we don't generate interior segment first, we need to
						//	ensure that vertex is exterior.

						for (vi=0;vi<ns;vi++)
							if (!used[vi] &&
								(wantInterior || !q.contains(vertices[fv[vi]]))) break;
						if (vi==ns) break;
						
						p0=vertices[fv[vi]];
						used[vi]=true;
						needPoint=false;
						onInter=false;
						};
			
					gf[0]=p0;				//	Store first point
					
					//	Loop, collecting all vertices of gon for this face
					
					collect:
					for (var tt5=0;tt5<10000;tt5++)
						{
						if (onInter)
							{
							//	We've just copied a point on the intersection of f1 and f2,
							//	specified by data index indx; move along the intersection
							//	line segment to the other point.
							
							eids[ngf++]=edgeFlag | f2;		//	Record edge ID
							indx=NDATA-indx;				//	Move to other point
							id=idF[f2];
							var phnum=id[indx+PH];
							if ((pn=ptsP[phnum][id[indx+PT]])==p0) break collect;
							gf[ngf]=pn;			//	Store new point
							
							if (phnum==ip)
								{
								//	We've hit an edge of our own polyhedron;
								//	we need to follow that edge.
							
								onInter=false;
								ie=id[indx+EDGE];
								ien=id[indx+INTER];
								vi=-1;
								needPoint=false;
								}
							else
								{
								//	Move to the other face-defined line segment for
								//	this intersection point.
							
								f2=id[indx+FACE];
								indx=id[indx+INDX];
								};
							}
						else if (vi>=0)
							{
							if (needPoint)
								{
								//	We just reached vertex vi of this face by
								//	travelling along the edge; we need to copy
								//	the vertex then see what follows.
								
								if ((pn=vertices[fv[vi]])==p0) break collect;
								gf[ngf]=pn;
								needPoint=false;
								used[vi]=true;
								};
						
							//	We just copied vertex # vi of current face

							eids[ngf++]=GonND.EDGE_ORIGINAL;	//	Original-edge ID
							
							//	If edge has no intersections, just go to next vertex.
							//	If edge has intersection(s), choose first one CCW around
							//	this face.

							ie=fe[vi];				//	Polyhedron edge number
							var ni=eIdata[ie][0];	//	Count of its intersections
							if (ni==0) vi=(vi+1)%ns;
							else
								{
								ien=(ni==1)?0:(p.face_edge_rev[f1][vi] ? 1 : 0);
								vi=-1;
								};
							needPoint=true;
							}
						else
							{
							if (needPoint)
								{
								//	We just reached intersection #ien of edge ie, by
								//	travelling along the edge; we need to copy the point
								//	and leave the edge.
							
								//	Get face for this edge intersection
								
								f2=eIdata[ie][2*ien+1];
								id=idF[f2];
								
								//	Ensure correct data index
								
								if (id[indx+PH]!=ip || id[indx+EDGE]!=ie) indx=NDATA-indx;
								
								//	Copy intersection point
								
								if ((pn=ptsP[ip][id[indx+PT]])==p0) break collect;
								gf[ngf]=pn;
								needPoint=false;
								
								//	Now pursue face-face intersection segment
								
								onInter=true;
								}
							else
								{
								//	We just copied intersection # ien of edge ie, after
								//	arriving on edge via face-face segment
							
								eids[ngf++]=GonND.EDGE_ORIGINAL;	//	Original-edge ID
								
								var nei=eIdata[ie][0];		//	Number of edge-intersections
								
								if (nei==2 && pass==0)
									{
									//	If this edge has two intersections, move to the
									//	other intersection immediately.

									ien=1-ien;
									}
								else
									{
									//	Identify vertex/edge number
								
									for (vi=0;vi<ns;vi++) if (fe[vi]==ie) break;
									
									//	If edge has only one intersection, or if we're
									//	on the later of two, go to the following vertex;
									//	otherwise go to the other intersection
									
									if (nei==1 || ien==(p.face_edge_rev[f1][vi]?0:1))
										vi=(vi+1)%ns;
									else ien=1-ien;
									};
								needPoint=true;
								};
							};
						};

					//	Copy polygon we accumulated in gf to output list
					
					var gnd=new Gon3D(gf,0,ngf);
					gnd.surfaceID=f1;
					gnd.flags=phFlag | pieceFlag | faceFlag;
					gnd.edgeIDs=new Array(ngf);
					for (var qq=0;qq<ngf;qq++) gnd.edgeIDs[qq]=eids[qq];
			 
					if (pass!=0 || wantInterior)
						{
						gons.push(gnd);
						};
					
					if (pass==0)
						{
						interiorPiece=gnd;		//	Make this available for later calcs
						break;					//	Single interior piece
						};
					};
				};
			}
		else
			{
			
			//	No face of other polyhedron cuts this face at all;
			//	so this face is either wholly inside or wholly outside
			
			var pieceFlag = q.contains(vertices[fv[0]]) ? Polyhedron.INTERIOR_PIECE : Polyhedron.EXTERIOR_PIECE;
			if ((pieceFlag & choiceFlags)!=0)
				{
				var gnd=p.faceToGon3D(f1);
				gnd.surfaceID=f1;
				gnd.flags=phFlag | pieceFlag | faceFlag;
				gons.push(gnd);
				};
			};
		};	//	for all faces
	};	//	for all polyhedra
return gons;
}

//	Return a specified face as a "Gon3D" object, with surfaceID equal to the
//	face number, and edgeIDs all indicating originals.

Polyhedron.prototype.faceToGon3D =
function (face)
{
var fv=this.face_vertices[face], ns=fv.length;
var g=new Array(ns);
for (var i=0;i<ns;i++) g[i]=this.vertices[fv[i]];
var gnd=new Gon3D(g,ns);

gnd.surfaceID=face;
gnd.setEdgeIDs();
return gnd;
}

//	Construct a convex polyhedron from a list of its face planes
//
//	Each face plane is described by an array of length 4, with the first
//	three components the outwards-directed unit normal to the face, and
//	the fourth component the value of the dot product with the normal
//	for the plane of the face.
//
//	faces			gives the face planes
//	nf					is the number of faces
//	toler				is the tolerance to use in determining if two points
//						are equal
//
//	Note that if any face planes fail to intersect the polyhedron, a new array
//	for the face data will be created with a smaller number of faces defined,
//	so the calling program should not assume that the original array is retained
//	by the Polyhedron object.

function PolyhedronFP(faces, nf, toler, makeEdgeArrays, makeVertexArrays)
{
this.faces=faces;
this.NF=nf;

this.face_vertices=new Array(nf);
this.face_edges=new Array(nf);
this.face_edge_rev=new Array(nf);

//	"bad" faces are those that are either parallel to other faces,
//	or those which don't intersect the convex hull of the other faces

var badFaceCount=0;
var badFaces=Utils.multiDim([nf],false);

//	Completed faces have had all their edges identified.

var completedFaceCount=0;
var completedFaces=new Array(nf);

//	Compute mutual dot products between face normals
//	fdp0 stores only dot product
//	fdp stores dot product tagged with face number, for sorting

var fdp0=Utils.multiDim([nf,nf],0), fdp=Utils.multiDim([nf,nf,2],0);
for (var i=0;i<nf;i++)
	{
	var fi=faces[i];
	for (var j=0;j<=i;j++)
		{
		var dp=0;
		for (var q=0;q<3;q++) dp+=fi[q]*faces[j][q];
		fdp0[i][j]=fdp0[j][i]=fdp[i][j][0]=fdp[j][i][0]=dp;
		fdp[i][j][1]=j;
		fdp[j][i][1]=i;
		
		//	Near-parallel normals are pathological, so we remove
		
		if (!badFaces[j] && j!=i && dp>1.0-toler) {badFaces[j]=true; badFaceCount++;};
		};
	};
	
//	edgeData holds data for edges defined by pairs of faces; entries are:
//
//		0	vertex number for start of edge
//		1	vertex number for end of edge
//		2	edge number
//		3	k such that fpV[i][j][k] refers to start of edge, and 1-k to end of
//			edge

var edgeData=Utils.multiDim([nf,nf],null);

//	lastEdge= 1+(other face number) for latest edge created for this face

var lastEdge=Utils.multiDim([nf],0);

//	nonEdgeFacePair is set true whenever we find that the line of intersection
//	of two faces definitely can't form an edge of the polyhedron

var nonEdgeFacePair=Utils.multiDim([nf,nf],0);

//	fpV holds data about vertices, indexed by face pairs.
//	The first two indices are face numbers, the second is 0 or 1 for
//	the two possible vertices shared by this face pair, then the first entry
//	is a vertex number, followed by a list of the other faces that intersect there.
//	Note that fpV[i][j] not being null doesn't mean face i and j have
//	an *edge* in common, as they might share only one vertex.

var fpV=Utils.multiDim([nf,nf],null);

//	faceToProcess gives face numbers in order in which we'll handle them
//	Initially it's set to 0,1,2,...nf, but we bring later entries ahead when
//	that allows us to keep processing faces where we already have at least one
//	edge

var faceToProcess=new Array(nf);

var vec=new Array(3), pt0=new Array(3), pt=null;
var ne=0, nv=0, startV=0, endV=0, firstV=0;
var vv=[];
var fv=Utils.multiDim([nf],0), fe=Utils.multiDim([nf],0), loF=Utils.multiDim([nf],0), hiF=Utils.multiDim([nf],0);

//	Candidates to look at for second plane intersecting current one
//	to make an edge

var candidates=Utils.multiDim([nf],0), nCand=0;

for (var l=0;l<nf;l++) faceToProcess[l]=l;

iLoop:
for (var ii=0;ii<nf;ii++)			//	For each face
	{
	var i=faceToProcess[ii];
	
	//	If there's any prospect of processing a face with some edge data already
	//	existing, try to do so
	
	if (lastEdge[i]==0 && ii<nf-1)
		{
		for (var l=ii+1;l<nf;l++)
			{
			var zz=faceToProcess[l];
			if (lastEdge[zz]>0)
				{
				faceToProcess[l]=i;		//	Swap current face to later in list
				i=zz;
				break;
				};
			};
		};
	
	if (badFaces[i]) continue;

	var fi=faces[i];			//	Normal vector and dot product value
									
	if (lastEdge[i]==0)
		{
		//	We have no record of another face that intersects face i.
		
		if (completedFaceCount>0)
			{
			//	We have at least one completed face, so if we can't find any more
			//	faces with existing edge data, that means all remaining faces must be
			//	bad ones that don't intersect the convex hull.
			
			badFaces[i]=true;
			badFaceCount++;
			for (var l=ii+1;l<nf;l++)
				{
				var zz=faceToProcess[l];
				if (!badFaces[zz])
					{
					badFaces[zz]=true;
					badFaceCount++;
					};
				};
			break iLoop;
			}
		else
			{
			//	We need to sort all faces by dot product, then look for intersections
			//	first between nearly (but non-)parallel faces
			
			Utils.qs(0,fdp[i],0,nf-1);
			nCand=0;
			for (var l=nf-2;l>=0;l--) candidates[nCand++]=Math.floor(fdp[i][l][1]);
			};
		}
	else
		{
		//	We have a record of another face that intersects face i, so
		//	make that our only candidate for creating the first edge.
		
		nCand=1;
		candidates[0]=lastEdge[i]-1;
		};

	var nfe=0;						//	Count of edges processed for this face
	startV=-1;						//	No previous vertex number known yet
	for (var tt6=0;tt6<10000;tt6++)					//	Loop until we've identified all edges
		{
		var dp=0, det=0;
		var j=0;
		
		var k=0;
		kLoop:
		for (k=0;k<nCand;k++)
			{
			j=candidates[k];
			if (edgeData[i][j]!=null) break kLoop;	//	Already have data for this edge
			
			if (badFaces[j] || nonEdgeFacePair[i][j])
				continue kLoop;			//	We know this i,j are no good
			
			dp=fdp0[i][j];
			det=1.0-dp*dp;
			if (det<toler)
				{
				//	Faces are parallel, so this i,j are no good
				
				nonEdgeFacePair[i][j]=nonEdgeFacePair[j][i]=true;
				continue kLoop;
				};
				
			//	Compute vector for line of intersection of
			//	face i and face j

			var fj=faces[j];
			Utils.cross2(fi,fj,vec);

			//	Chop the line of intersection down to a line segment
			//	by imposing inequalities from the requirement that the edge lies on
			//	one side of every face plane (other than i and j, and bad faces to be
			//	discarded)
			
			var loPar=Number.NEGATIVE_INFINITY;
			var hiPar=Number.POSITIVE_INFINITY;
			var nLo=0, nHi=0;
			endV=-1;
			
			//	If we have pre-existing vertex info, make use of that.
			
			var chopFirst=Utils.multiDim([nf],false);
			var mv=fpV[i][j], mvSign=1;
			if (mv!=null)
				{
				//	If edge direction agrees with one of the outwards face normals,
				//	vertex lies at the *end* of the edge, rather than the beginning
				
				if (Utils.dot(vec,faces[mv[0][1]])>0) mvSign=-1;
				
				for (var l=0;l<2;l++)	//	Loop for up to 2 different entries
					{
					var mvl=mv[l];
					if (mvl!=null)
						{
						var pp=mvl.length;
						if (l==0 && mvSign==1 || l==1 && mvSign==-1)
							{
							if (startV>0 && startV!=mvl[0]) continue kLoop;
							startV=mvl[0];
							for (var m=1;m<pp;m++)
								chopFirst[loF[nLo++]=mvl[m]]=true;
							}
						else
							{
							endV=mvl[0];
							for (var m=1;m<pp;m++)
								chopFirst[hiF[nHi++]=mvl[m]]=true;
							};
						};
					};
				if (endV==startV) continue kLoop;
				};
			
			var startUnchopped=startV<0;
			var endUnchopped=endV<0;
			
			if (startUnchopped || endUnchopped)
				{
				if (startUnchopped)
					{
					//	Compute point pt on line of
					//	intersection of face i and face j
					
					var c1=fi[3], c2=fj[3];
					var a=(c1-dp*c2)/det, b=(c2-dp*c1)/det;
					for (var l=0;l<3;l++) pt0[l]=a*fi[l]+b*fj[l];
					pt=pt0;
					}
				else
					{
					pt=vv[startV];
					loPar=0.0;
					};

				//	Now chop the line
				
				var chopFaces=Utils.multiDim([nf],false);
				var nc=0;
				for (var l=0;l<nf;l++)
					if (chopFirst[l]) chopFaces[nc++]=l;
				for (var l=0;l<nf;l++)
					if (!chopFirst[l] && !badFaces[l] && l!=i && l!=j && !completedFaces[l])
						chopFaces[nc++]=l;

				for (var ic=0;ic<nc;ic++)
					{
					var f=chopFaces[ic];
					var ff=faces[f];
					var dpt=Utils.dot(ff,pt), dvec=Utils.dot(ff,vec);
					if (Math.abs(dvec)<toler)
						{
						//	Plane is parallel to line of intersection, so it either
						//	rules out whole line, or has no effect
						
						if (dpt<=ff[3]+toler) continue;		//	Line is on right side of plane
						else
							{
							//	Line is on wrong side of plane, so this i,j are no good
							nonEdgeFacePair[i][j]=nonEdgeFacePair[j][i]=true;
							continue kLoop;
							}
						}
						
					//	Work out parameter where we hit plane, and adjust
					//	loPar or hiPar accordingly to chop to smaller segment
					
					var par=(ff[3]-dpt)/dvec;
					if (dvec>0)
						{
						if (par<hiPar-toler)
							{
							hiPar=par;
							if (endUnchopped) {hiF[0]=f; nHi=1;};
							}
						else if (par<hiPar+toler && endUnchopped) hiF[nHi++]=f;
						}
					else
						{
						if (par>loPar+toler)
							{
							loPar=par;
							if (startUnchopped) {loF[0]=f; nLo=1;};
							}
						else if (par>loPar-toler && startUnchopped) loF[nLo++]=f;
						};
					
					//	Give up on this i,j if we've chopped edge down to nothing
					
					if (hiPar<loPar+toler)
						{
						nonEdgeFacePair[i][j]=nonEdgeFacePair[j][i]=true;
						continue kLoop;
						};
					};	//	End loop for chopping faces
				};
				
			//	If edge is unchopped, we have pathological, unclosed polyhedron
			
			if (nLo==0 || nHi==0)
				{
				this.NF=0;
				return;
				};
				
			//	Edge has survived; now record some details
			
			var eij=edgeData[i][j]=[0,0,0,0];
			var eji=edgeData[j][i]=[0,0,0,0];

			eij[2]=eji[2]=ne++;
			
			//	Record vertex data for start of edge
			
			if (startUnchopped)
				{
				var sv=new Array(3);
				for (var l=0;l<3;l++) sv[l]=pt[l]+loPar*vec[l];
				startV=nv++;
				vv.push(sv);
				loF[nLo++]=i;
				loF[nLo++]=j;
				for (var l=0;l<nLo-1;l++)
				for (var m=l+1;m<nLo;m++)
					{
					var mvd=Utils.multiDim([nLo-1],0);
					mvd[0]=startV;
					var ll=loF[l], mm=loF[m], nn=0;
					if (!fpV[ll][mm]) fpV[ll][mm]=fpV[mm][ll]=[null,null];
					else nn=1;
					fpV[ll][mm][nn]=mvd;
					var zz=1;
					for (var n=0;n<nLo;n++)
						if (n!=l && n!=m) mvd[zz++]=loF[n];
					};
				};
			eij[0]=eji[1]=startV;

			//	Record vertex data for end of edge
			
			if (endUnchopped)
				{
				var ev=new Array(3);
				for (var l=0;l<3;l++) ev[l]=pt[l]+hiPar*vec[l];
				endV=nv++;
				vv.push(ev);
				hiF[nHi++]=i;
				hiF[nHi++]=j;
				for (var l=0;l<nHi-1;l++)
				for (var m=l+1;m<nHi;m++)
					{
					var mvd=Utils.multiDim([nHi-1],0);
					mvd[0]=endV;
					var ll=hiF[l], mm=hiF[m], nn=0;
					if (!fpV[ll][mm]) fpV[ll][mm]=fpV[mm][ll]=[null,null];
					else nn=1;
					fpV[ll][mm][nn]=mvd;
					var zz=1;
					for (var n=0;n<nHi;n++)
						if (n!=l && n!=m) mvd[zz++]=hiF[n];
					};
				};
			eij[1]=eji[0]=endV;
			
			eji[3]=1-(
				eij[3]=(fpV[i][j][0][0]==startV) ? 0 : 1
				);

			break kLoop;
			};	//	End kLoop, searching for second face

		if (k==nCand || nfe==nf)
			{
			if (!badFaces[i])
				{
				badFaces[i]=true;
				badFaceCount++;
				};
			continue iLoop;
			};
			
		var eij=edgeData[i][j];

		//	Remember vertex at start of first edge
		
		if (nfe==0) firstV=eij[0];
		
		//	Note face in lastEdge so we can start face j here, rather than
		//	searching
		
		lastEdge[j]=1+i;
		
		//	Log vertex and edge numbers for this edge of the face
		
		fv[nfe]=eij[0];
		fe[nfe]=eij[2];
		
		nfe++;					//	Count up edges created for current face
	
		//	Remember end vertex number of this edge
		//	to use as start vertex number of next edge
		
		startV=eij[1];
		if (startV==firstV) break;	//	We've closed the face, so we're done
		
		//	Make faces that cut end of edge candidates for next face to intersect
		//	current face, face i
		
		var fl=fpV[i][j][1-eij[3]];
		nCand=fl.length-1;
		for (var l=0;l<nCand;l++) candidates[l]=fl[l+1];
		};	//	End loop finding edges for face i
		
	completedFaceCount++;
	completedFaces[i]=true;
	
	var fvi=this.face_vertices[i]=new Array(nfe);
	var fei=this.face_edges[i]=new Array(nfe);
	for (var l=0;l<nfe;l++)
		{
		fvi[l]=fv[l];
		fei[l]=fe[l];
		};
		
	};	//	End i-loop scanning all faces

//	We encountered faces that don't intersect the polyhedron, so they have
//	to be removed from the list
	
if (badFaceCount>0)
	{
	var nfn=nf-badFaceCount;
	var tmp1=new Array(nfn);
	var tmp2=new Array(nfn);
	var tmp3=new Array(nfn);
	var i=0;
	for (var j=0;j<nf;j++)
		if (!badFaces[j])
			{
			tmp1[i]=faces[j];
			tmp2[i]=this.face_vertices[j];
			tmp3[i]=this.face_edges[j];
			i++;
			};
	this.faces=faces=tmp1;
	this.NF=nf=i;
	this.face_vertices=tmp2;
	this.face_edges=tmp3;
	};
	
//	Create arrays of edge data

this.NE=ne;
if (makeEdgeArrays)
	{
	this.edge_vertices=new Array(ne);
	this.edge_faces=Utils.multiDim([ne,2],0);
	for (var i=0;i<nf;i++)
		{
		var nfe=this.face_edges[i].length;
		this.face_edge_rev[i]=new Array(nfe);
		for (var j=0;j<nfe;j++)
			{
			var en=this.face_edges[i][j];
			var v1=this.face_vertices[i][j], v2=this.face_vertices[i][(j+1)%nfe];
			if (this.edge_vertices[en]==null)
				{
				var evn=this.edge_vertices[en]=new Array(2);
				evn[0]=v1;
				evn[1]=v2;
				this.edge_faces[en][0]=i;
				}
			else
				{
				this.face_edge_rev[i][j]=true;
				this.edge_faces[en][1]=i;
				};
			};
		};
	};
//	Create array of vertices and vertex data

this.NV=nv;
this.vertices=Utils.multiDim([nv,3],0);
for (var i=0;i<nv;i++) this.vertices[i]=vv[i];

if (makeVertexArrays)
	{
	this.vertex_vertices=new Array(nv);		//	List of vertices around each vertex
	this.vertex_edges=new Array(nv);			//	List of edges around each vertex
	this.vertex_faces=new Array(nv);			//	List of faces around each vertex
	};

this.isometries=null;
this.even_isometry=null;
this.vertex_permutations=null;
this.face_permutations=null;
}

PolyhedronFP.prototype = Polyhedron.prototype;
    