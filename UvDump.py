import bpy
import decimal
from bpy_extras.io_utils import ExportHelper

class UvDumpWriteData(bpy.types.Operator, ExportHelper):
    bl_idname="uvdump.write_data"
    bl_label="Write UV data to file"
    filename_ext=""
    
    def execute(self, context):
        print(self.filepath)
        uv_dump(self.filepath)
        return {'FINISHED'}

# create a new context for this task
ctx = decimal.Context()

# 20 digits should be enough for everyone :D
ctx.prec = 20

def float_to_str(f):
    """
    Convert the given float to a string,
    without resorting to scientific notation
    """
    d1 = ctx.create_decimal(repr(f))
    return format(d1, 'f')

def uv_dump(filename):
    obj = bpy.context.active_object.data
    obj_mat = bpy.context.active_object.matrix_world
    uv_layer = obj.uv_layers.active.data

    vertices = [None] * len(obj.vertices)
    polygons = [[] for j in range(0, len(obj.polygons))]
    print(polygons)

    for v in range(0, len(obj.vertices)):
        vert_vec = obj_mat @ obj.vertices[v].co
        vertices[v] = (vert_vec[0], vert_vec[1], vert_vec[2])

    for poly in obj.polygons:
        print("Polygon index: %d, length: %d" % (poly.index, poly.loop_total))
        for loop_index in range(poly.loop_start, poly.loop_start + poly.loop_total):
            print("    Vertex: %d" % obj.loops[loop_index].vertex_index)
            print("    UV: %r" % uv_layer[loop_index].uv)
            
            vertex_idx = obj.loops[loop_index].vertex_index
            uv_vec = uv_layer[loop_index].uv

            polygons[poly.index].append((obj.loops[loop_index].vertex_index, (uv_vec[0], uv_vec[1])))

    for poly in polygons:
        print(poly)
    print(vertices)
        
    curves = list(filter(lambda o: o.type == "CURVE" and o != bpy.context.active_object, bpy.context.selected_objects))
    splines = []
    for curve in curves:
        for spline in curve.data.splines:
            spline_points = []
            for bezier_point in spline.bezier_points:
                _left = curve.matrix_world @ bezier_point.handle_left
                _center = curve.matrix_world @ bezier_point.co
                _right = curve.matrix_world @ bezier_point.handle_right
                handle_left = (_left[0], _left[1], _left[2])
                loc = (_center[0], _center[1], _center[2])
                handle_right = (_right[0], _right[1], _right[2])
                spline_points.append((handle_left, loc, handle_right))
            splines.append(spline_points)
                
    with open(filename, "w") as file_out:
        file_out.write(str(len(vertices)) + '\n')
        for i in range(0, len(vertices)):
            file_out.write(float_to_str(vertices[i][0]) + ' ' + float_to_str(vertices[i][1]) + ' ' + float_to_str(vertices[i][2]) + '\n')
        file_out.write(str(len(polygons)) + '\n')
        for i in range(0, len(polygons)):
            for j in range(0, len(polygons[i])):
                file_out.write(str(polygons[i][j][0]) + ' ');
                file_out.write('(' + float_to_str(polygons[i][j][1][0]) + ', ' + float_to_str(polygons[i][j][1][1]) + ') ')
            file_out.write('\n')
        file_out.write(str(len(splines)) + '\n')
        for spline in splines:
            for point in spline:
                tuple_format = lambda t: "({}, {}, {})".format(float_to_str(t[0]), float_to_str(t[1]), float_to_str(t[2]))
                file_out.write("({}, {}, {}) ".format(tuple_format(point[0]), tuple_format(point[1]), tuple_format(point[2])))
            file_out.write('\n')
            
bpy.utils.register_class(UvDumpWriteData)
bpy.ops.uvdump.write_data('INVOKE_DEFAULT', filepath="uv_data")
