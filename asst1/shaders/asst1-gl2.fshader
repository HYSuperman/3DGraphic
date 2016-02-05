uniform float uVertexScale;
uniform float uTime;
uniform sampler2D uTexUnit0, uTexUnit1;

varying vec2 vTexCoord0, vTexCoord1;
varying vec3 vColor;

void main(void) {
  vec4 color = vec4(vColor.x, vColor.y, vColor.z, 1);
  vec4 texColor0 = texture2D(uTexUnit0, vTexCoord0);
  vec4 texColor1 = texture2D(uTexUnit1, vTexCoord1);

  float lerper = clamp(.5 *uVertexScale, 0., 1.);
  float lerper2 = clamp(.5 * uVertexScale + 1.0, 0.0, 1.0);
  
  float timme = sin(uTime);
  if(timme < 0.0)
    lerper = clamp(timme, 0.0, 1.0);
  else
    lerper = clamp(-timme, 0.0, 1.0);
  lerper2 = 1.0;
/*
  lerper = clamp(timme, 0., 1.);
  if(timme > 0.0)
  lerper2 = clamp(timme + 1.0, 0.0, 1.0);
  else
  lerper2 = clamp(-timme, 0.0, 1.0);
  */
  // ======================================
  // TODO: use sin and uTime to interpolate
  //       between the two images
  // ======================================

  //gl_FragColor = ((lerper)*texColor1 + (1.0-lerper)*texColor0) * lerper2 + color * (1.0-lerper2);
  gl_FragColor = ((lerper)*texColor1 + (1.0-lerper)*texColor0)*lerper2 ;
}
