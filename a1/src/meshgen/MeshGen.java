package meshgen;

import java.util.ArrayList;
import java.util.Arrays;
import java.io.IOException;
import math.Vector2;
import math.Vector3;
import meshgen.OBJMesh;
import meshgen.OBJFace;

public class MeshGen{

	public static void main (String[] args) throws IOException{
		//default values of n,m,r
		int n = 32;
		int m = 16;
		double r = 0.25;
		String fileName = "";
		String outFile = "";
		
		if (args[0].equals("-g")){
			if(args[1].equals("sphere")){
				if(args.length >= 3){
					if(args[2].equals("-n")){
						n = Integer.parseInt(args[3]);
						if(args[4].equals("-m")){
							m = Integer.parseInt(args[5]);
							fileName = args[7];
						}
						else if (args[4].equals("-o")){fileName=args[5];}
					}
					else if(args[2].equals("-m")){
						m = Integer.parseInt(args[3]);
						fileName = args[5];
					}
					else if(args[2].equals("-o")){
						fileName = args[3];
					}
				}
				createSphere(n,m,fileName);
			}
			else if(args[1].equals("cylinder")){
				if(args[2].equals("-n")){
					n = Integer.parseInt(args[3]);
					fileName = args[5];
				}
				else if(args[2].equals("-o")){
					fileName = args[3];
				}
				createCylinder(n,fileName);
			}
			else if(args[1].equals("torus")){
				if(args.length >= 3){
					if(args[2].equals("-r")){
						r = Double.parseDouble(args[3]);
						fileName = args[5];
					}
					else{fileName = args[3];}
				}
				// createTorus(r, fileName);
			}
			else{System.out.println("Choose cylinder, sphere or torus");}
		}
		else{
			fileName = args[1];
			outFile = args[3];
			createNormals(fileName, outFile);
		}
		
		// createNormals("horse-nonorms.obj", "test.obj");
	}

	public static void createCylinder(int n, String f) throws IOException{
		double innerAngle = (2*Math.PI)/n;
		OBJMesh cyl = new OBJMesh();
		ArrayList<Vector3> cylinderVertex = cyl.positions;
		ArrayList<Vector3> cylinderNormals = cyl.normals;
		ArrayList<Vector2> cylinderTexture = cyl.uvs;
		ArrayList<OBJFace> cylinderFaces = cyl.faces;

		//creates positions and normals
		for(int i = 0; i < 2*n; i++){
			float x_coord, z_coord;
			if(i < n){
				x_coord = (float) Math.cos(innerAngle*i);
				z_coord = (float) Math.sin(-innerAngle*i);
				cylinderVertex.add(new Vector3(x_coord, (float)1.0, z_coord));
				cylinderNormals.add(new Vector3(x_coord, 0, z_coord));
			}
			else{
				x_coord = (float) Math.cos(innerAngle*(i-n));
				z_coord = (float) Math.sin(-innerAngle*(i-n));
				cylinderVertex.add(new Vector3(x_coord, (float)-1.0, z_coord));
			}
		}
		cylinderVertex.add(new Vector3((float)0.0, (float)1.0, (float)0.0));
		cylinderVertex.add(new Vector3((float)0.0, (float)-1.0, (float)0.0));
		cylinderNormals.add(new Vector3((float)0.0, (float)1.0,(float)0.0));
		cylinderNormals.add(new Vector3((float)0.0, (float)-1.0, (float)0.0));
		//creates texture coordinates
		//int c1 = 0;
		//int c2 = 0;
		//int c3 = 0;
		for(int i = 0; i <= n; i++){
			cylinderTexture.add(new Vector2((float)i/n, (float)0.5));
			//c1++;
		}
		for(int i = 0; i <= n;i++){
				cylinderTexture.add(new Vector2((float)i/n, (float)0.0));
				//c1++;
		}
		
		for(int i = 0; i < n; i++) {
				//texture coordinates for bottom
				cylinderTexture.add(new Vector2((float)(0.25+0.25*Math.cos(innerAngle*i)),
						(float)(0.75+0.25*Math.sin(innerAngle*i))));
				//c2++;
			}
		
		for(int i = 0; i < n; i++) {
				//texture coordinates for top
				cylinderTexture.add(new Vector2((float)(0.75+0.25*Math.cos(innerAngle*i)), 
						(float)(0.75+0.25*Math.sin(innerAngle*i))));
				//c3++;
			}
		
		//System.out.println("There are "+c1+" coordinates for the sides, "+c2+" coordinates for the bottom, and "+c3+" coordinates for the top");
		cylinderTexture.add(new Vector2((float)0.75,(float)0.75));
		cylinderTexture.add(new Vector2((float)0.25,(float)0.75));
		//might have to increase all indicies by 1
		for(int i = 0; i < n-1; i++){
			OBJFace topTriangle = new OBJFace(3, true, true);
			OBJFace side1 = new OBJFace(3, true, true);
			OBJFace side2 = new OBJFace(3, true, true);
			OBJFace bottomTriangle = new OBJFace(3, true, true);
			
			topTriangle.setVertex(0, i, 3*n+i+2, n);
			topTriangle.setVertex(1, i+1, 3*n+i+3, n);
			topTriangle.setVertex(2, 2*n, 4*n+2, n);
			cylinderFaces.add(topTriangle);
			
			side1.setVertex(0, i, i, i);
			side1.setVertex(1, i+n, i+n+1, i);
			side1.setVertex(2, i+1, i+1, i+1);
			cylinderFaces.add(side1);
			
			side2.setVertex(0, i+1, i+1, i+1);
			side2.setVertex(1, i+n, i+n+1, i);
			side2.setVertex(2, i+n+1, i+n+2, i+1);
			cylinderFaces.add(side2);
	
			//MUST FIX
			bottomTriangle.setVertex(0, i+n, 2*n+i+2, n+1);
			bottomTriangle.setVertex(1, 2*n+1, 4*n+3, n+1);
			bottomTriangle.setVertex(2, i+n+1, 2*n+i+3, n+1);
			cylinderFaces.add(bottomTriangle);
		}
		OBJFace lastTop = new OBJFace(3, true, true);
		lastTop.setVertex(0, n-1, 4*n, n);
		lastTop.setVertex(1, 0, 4*n+1, n);
		lastTop.setVertex(2, 2*n, 4*n+2, n);
		cylinderFaces.add(lastTop);
		
		OBJFace lastSide1 = new OBJFace(3, true, true);
		lastSide1.setVertex(0, n-1, n-1, n-1);
		lastSide1.setVertex(1, 2*n-1, 2*n, n-1);
		lastSide1.setVertex(2, 0, n, 0);
		cylinderFaces.add(lastSide1);
		
		OBJFace lastSide2 = new OBJFace(3, true, true);
		lastSide2.setVertex(0, 0, n, 0);
		lastSide2.setVertex(1, 2*n-1, 2*n, n-1);
		lastSide2.setVertex(2, n, 2*n+1, 0);
		cylinderFaces.add(lastSide2);
		
		OBJFace lastBottom = new OBJFace(3, true, true);
		lastBottom.setVertex(0, n, 3*n+1, n+1);
		lastBottom.setVertex(1, 2*n-1, 3*n, n+1);
		lastBottom.setVertex(2, 2*n+1, 4*n+3, n+1);
		cylinderFaces.add(lastBottom);

		cyl.isValid(true);
		cyl.writeOBJ(f);
	}

	public static void createSphere(int u, int v, String f) throws IOException{
		OBJMesh sph = new OBJMesh();
		ArrayList<Vector3> sphereVertex = sph.positions;
		ArrayList<Vector3> sphereNormal = sph.normals;
		ArrayList<Vector2> sphereTexture = sph.uvs;
		ArrayList<OBJFace> sphereFaces = sph.faces;
		sphereTexture.add(new Vector2((float) 0, (float) 0));
		double sideAngle = 2*Math.PI/u;
		double vertAngle = Math.PI/v;
		for(int i = 1; i < v; i++){
			float x_coord,y_coord, z_coord, r;
			y_coord = (float) Math.sin(vertAngle*i+ Math.PI/2);
			r = (float) Math.cos(vertAngle*i+Math.PI/2);
			for(int j = 0; j < u; j++){
				x_coord = r * (float) Math.cos(sideAngle * j);
				z_coord = r * (float) Math.sin(sideAngle * j);
				sphereVertex.add(new Vector3(x_coord, y_coord, z_coord));
				sphereNormal.add(new Vector3(x_coord, y_coord, z_coord));
			}
		}
		
		//top pole
		sphereVertex.add(new Vector3((float)0.0, (float)1.0, (float)0.0));
		sphereNormal.add(new Vector3((float)0.0, (float)1.0, (float)0.0));
		
		//bottom pole
		sphereVertex.add(new Vector3((float)0.0, (float)-1.0, (float)0.0));
		sphereNormal.add(new Vector3((float)0.0, (float)-1.0, (float)0.0));
		
		//create texture verticies
		
		// pole #1 column #1
		for(int i = 1; i < v + 1; i++) {
			sphereTexture.add(new Vector2 (0, ((float)i)/v));
		}
		
		// center part
		for (int i = 0; i < v + 1; i ++) {
			for (int j = 1; j < u; j ++)
				sphereTexture.add(new Vector2 (((float)i)/u, ((float)j)/v));
		}
		
		for (int i = 0; i < v; i ++) {
			
		}
		
		//create faces
		for(int i = 0; i < u-1; i++){
			for(int j = 0; j < v-2; j++){
				OBJFace t1 = new OBJFace(3,true,true);
				t1.setVertex(2, i+j*u, i+j*u, i+j*u);
				t1.setVertex(1, i+(j+1)*u, i+(j+1)*u, i+(j+1)*u);
				t1.setVertex(0, i+j*u+1, i+j*u+1, i+j*u+1);
				sphereFaces.add(t1);
				OBJFace t2 = new OBJFace(3,true,true);
				t2.setVertex(2, i+j*u+1, i+j*u+1, i+j*u+1);
				t2.setVertex(1, i+(j+1)*u, i+(j+1)*u, i+(j+1)*u);
				t2.setVertex(0, i+(j+1)*u+1, i+(j+1)*u+1, i+(j+1)*u+1);
				sphereFaces.add(t2);
			}
		}
		//close surface
		for(int i = 1; i < v-1; i++){
			OBJFace f1 = new OBJFace(3,true,true);
			//texture coordinates are different here
			f1.setVertex(0, i*u, i*u, i*u);
			f1.setVertex(1, i*u-1, i*u-1, i*u-1);
			f1.setVertex(2, (i-1)*u, (i-1)*u, (i-1)*u);
			sphereFaces.add(f1);

			OBJFace f2 = new OBJFace(3,true,true);
			//texture coordinates are different here
			f2.setVertex(2, i*u, i*u, i*u);
			f2.setVertex(1, i*u-1, i*u-1, i*u-1);
			f2.setVertex(0, (i+1)*u-1, (i+1)*u-1, (i+1)*u-1);
			sphereFaces.add(f2);
		}
		for(int i = 0; i < u-1; i++){
			OBJFace topTriag = new OBJFace(3,true, true);
			topTriag.setVertex(0, i, i, i);
			topTriag.setVertex(1, u*(v-1), u*(v-1), u*(v-1));
			topTriag.setVertex(2, i+1, i+1, i+1);
			sphereFaces.add(topTriag);
			
			OBJFace botTriag = new OBJFace(3,true,true);
			botTriag.setVertex(2, u*(v-2)+i, u*(v-2)+i, u*(v-2)+i);
			botTriag.setVertex(1, sphereVertex.size()-1, sphereVertex.size()-1, sphereVertex.size()-1);
			botTriag.setVertex(0, u*(v-2)+i+1, u*(v-2)+i+1, u*(v-2)+i+1);
			sphereFaces.add(botTriag);
		}
		OBJFace closeTop = new OBJFace(3,true, true);
		closeTop.setVertex(2, 0, 2, 0);
		closeTop.setVertex(1, u*(v-1), u*(v-1), u*(v-1));
		closeTop.setVertex(0, u-1, u-1, u-1);
		sphereFaces.add(closeTop);
		OBJFace closeBot = new OBJFace(3,true,true);
		closeBot.setVertex(0, (v-2)*u, (v-2)*u, (v-2)*u);
		closeBot.setVertex(1, (v-1)*u+1, (v-1)*u+1, (v-1)*u+1);
		closeBot.setVertex(2, (v-1)*u-1, (v-1)*u-1, (v-1)*u-1);		
		sphereFaces.add(closeBot);

		sph.isValid(true);
		sph.writeOBJ(f);
	}
 
	public static void createTor(int u, int v, String f) throws IOException{
		OBJMesh tor = new OBJMesh();
		ArrayList<Vector3> torVertex = tor.positions;
		ArrayList<Vector3> torNormal = tor.normals;
		ArrayList<Vector2> torTexture = tor.uvs;
		ArrayList<OBJFace> torFaces = tor.faces;
		double sideAngle = 2*Math.PI/u;
		double vertAngle = Math.PI/v;
		for(int i = 1; i < v; i++){
			float x_coord,y_coord, z_coord, r;
			y_coord = (float) Math.sin(vertAngle*i+ Math.PI/2);
			r = (float) Math.cos(vertAngle*i+Math.PI/2);
			for(int j = 0; j < u; j++){
				x_coord = (float) Math.cos(sideAngle * j);
				z_coord = (float) Math.sin(sideAngle * j);
				torVertex.add(new Vector3(x_coord, y_coord, z_coord));
				torNormal.add(new Vector3(x_coord, y_coord, z_coord));
			}
		}
		
		//top pole
		torVertex.add(new Vector3((float)0.0, (float)1.0, (float)0.0));
		torNormal.add(new Vector3((float)0.0, (float)1.0, (float)0.0));
		
		//bottom pole
		torVertex.add(new Vector3((float)0.0, (float)-1.0, (float)0.0));
		torNormal.add(new Vector3((float)0.0, (float)-1.0, (float)0.0));
		
		//create texture verticies
		
		// pole #1 column #1
		for(int i = 1; i < v + 1; i++) {
			torTexture.add(new Vector2 (0, ((float)i)/v));
		}
		
		// center part
		for (int i = 0; i < v + 1; i ++) {
			for (int j = 1; j < u; j ++)
				torTexture.add(new Vector2 (((float)i)/u, ((float)j)/v));
		}
		
		for (int i = 0; i < v; i ++) {
			
		}
		
		//create faces
		for(int i = 0; i < u-1; i++){
			for(int j = 0; j < v-2; j++){
				OBJFace t1 = new OBJFace(3,true,true);
				t1.setVertex(2, i+j*u, i+j*u, i+j*u);
				t1.setVertex(1, i+(j+1)*u, i+(j+1)*u, i+(j+1)*u);
				t1.setVertex(0, i+j*u+1, i+j*u+1, i+j*u+1);
				torFaces.add(t1);
				OBJFace t2 = new OBJFace(3,true,true);
				t2.setVertex(2, i+j*u+1, i+j*u+1, i+j*u+1);
				t2.setVertex(1, i+(j+1)*u, i+(j+1)*u, i+(j+1)*u);
				t2.setVertex(0, i+(j+1)*u+1, i+(j+1)*u+1, i+(j+1)*u+1);
				torFaces.add(t2);
			}
		}
		//close surface
		for(int i = 1; i < v-1; i++){
			OBJFace f1 = new OBJFace(3,true,true);
			//texture coordinates are different here
			f1.setVertex(0, i*u, i*u, i*u);
			f1.setVertex(1, i*u-1, i*u-1, i*u-1);
			f1.setVertex(2, (i-1)*u, (i-1)*u, (i-1)*u);
			torFaces.add(f1);

			OBJFace f2 = new OBJFace(3,true,true);
			//texture coordinates are different here
			f2.setVertex(2, i*u, i*u, i*u);
			f2.setVertex(1, i*u-1, i*u-1, i*u-1);
			f2.setVertex(0, (i+1)*u-1, (i+1)*u-1, (i+1)*u-1);
			torFaces.add(f2);
		}
		for(int i = 0; i < u-1; i++){
			OBJFace topTriag = new OBJFace(3,true, true);
			topTriag.setVertex(0, i, i, i);
			topTriag.setVertex(1, u*(v-1), u*(v-1), u*(v-1));
			topTriag.setVertex(2, i+1, i+1, i+1);
			torFaces.add(topTriag);
			
			OBJFace botTriag = new OBJFace(3,true,true);
			botTriag.setVertex(2, u*(v-2)+i, u*(v-2)+i, u*(v-2)+i);
			botTriag.setVertex(1, torVertex.size()-1, torVertex.size()-1, torVertex.size()-1);
			botTriag.setVertex(0, u*(v-2)+i+1, u*(v-2)+i+1, u*(v-2)+i+1);
			torFaces.add(botTriag);
		}
		OBJFace closeTop = new OBJFace(3,true, true);
		closeTop.setVertex(2, 0, 2, 0);
		closeTop.setVertex(1, u*(v-1), u*(v-1), u*(v-1));
		closeTop.setVertex(0, u-1, u-1, u-1);
		torFaces.add(closeTop);
		OBJFace closeBot = new OBJFace(3,true,true);
		closeBot.setVertex(0, (v-2)*u, (v-2)*u, (v-2)*u);
		closeBot.setVertex(1, (v-1)*u+1, (v-1)*u+1, (v-1)*u+1);
		closeBot.setVertex(2, (v-1)*u-1, (v-1)*u-1, (v-1)*u-1);		
		torFaces.add(closeBot);

		tor.isValid(true);
		tor.writeOBJ(f);
	}

	public static void createNormals(String infile, String outfile) throws IOException{
		//creates OBJMesh object out of input file
		OBJMesh inMesh = new OBJMesh();
		inMesh.parseOBJ(infile);
		//references to inMesh's positions, normals, faces arraylists
		ArrayList<Vector3> inPos = inMesh.positions;
		ArrayList<Vector3> inNorm = inMesh.normals; //will be empty
		ArrayList<OBJFace> inFace = inMesh.faces;
		//inserts empty Vector3 into each space for 
		for(int i = 0; i < inPos.size(); i++){
			inNorm.add(i, new Vector3((float)0.0, (float)0.0, (float)0.0));
		}
		System.out.println(inPos.size());
		System.out.println(inNorm.size());
		//takes each face
		for(int i = 0; i < inFace.size(); i++){
			OBJFace tri = inFace.get(i);
			tri.normals = tri.positions;
			//positions for the three points 
			Vector3 p1 = inPos.get(tri.positions[0]);
			Vector3 p2 = inPos.get(tri.positions[1]);
			Vector3 p3 = inPos.get(tri.positions[2]);
			//gets normal
			Vector3 triNorm = p1.clone().sub(p2).cross(p3.clone().sub(p2)).normalize();
			triNorm.negate();
			System.out.println("face " + i + ": normal is " + triNorm);
			for(int p = 0; p < 3; p++){
				if(tri.normals[p] < inPos.size()){
					inNorm.get(tri.normals[p]).add(triNorm);
				}
			}
		}
		for(int i = 0; i < inNorm.size(); i++){
			inNorm.get(i).normalize();
		}
		System.out.println(inNorm.size());
		//creates output of OBJMesh object that goes in output file
		if(inMesh.isValid(true)){
			inMesh.writeOBJ(outfile);
		}
	}
}