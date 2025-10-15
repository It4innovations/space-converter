#####################################################################################################################
# Copyright(C) 2023-2025 IT4Innovations National Supercomputing Center, VSB - Technical University of Ostrava
#
# This program is free software : you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#####################################################################################################################


from . import bspace_pref
from . import bspace_remote

import numpy as np
import bpy
import pathlib
from ctypes import *
import struct
import mathutils

import bpy
import gpu
from gpu_extras.batch import batch_for_shader
import blf

from bpy_extras import view3d_utils

##############################################################
# value_convert_items = [    
#     ("0", "DEFAULT", ""),
#     #("1", "LOG(V)", ""),    
#     #("2", "LOG(V/KPC^3)", "")
# ]

anim_items = [
    ("0", "NONE", ""),
    ("1", "ALL(path)", ""),
    ("2", "ALL(merge)", ""),
    ("3", "FRAME(extract)", ""),
    #("4", "FRAME(cache)", ""),
]

extracted_type_items = [    
    ("0", "Sparse", ""),
    ("1", "Dense", ""),
    ("2", "Particle", ""),
]

dense_type_items = [    
    #("0", "NONE", ""),
    ("1", "Cubic", ""),
    ("2", "Quintic", ""),
    ("3", "WendlandC2", ""),
    ("4", "WendlandC4", ""),
    ("5", "WendlandC6", ""),
    ("6", "WendlandC8", ""),
]

dense_norm_items = [    
    ("0", "NONE", ""),
    ("1", "Count", ""),
    ("2", "SPHInterpolation", ""),
]

file_type_items = [
    ("0", "NONE", ""),
    ("1", "OPENVDB", ""),
    ("2", "NANOVDB", ""),
    ("3", "PATH", ""),
    ("4", "RAW_PART", ""),
]

slice_axis_items = [
    ("x", "X", ""),
    ("y", "Y", ""),
    ("z", "Z", ""),
]
# ##############################################################
# Callback function to update the scale based on the bbox_size
def update_bbox_size(self, context):
    context.scene.bspace.set_bbox_size(context, self.bbox_size)

class BSPACE_PG_SETTINGS(bpy.types.PropertyGroup):    
    # local_temp_dir_path: bpy.props.StringProperty(
    #     name="Temp", 
    #     subtype='DIR_PATH', 
    #     default=''
    #     )  # type: ignore
    
    grid_dim: bpy.props.IntProperty(
        name="Grid Dim", 
        default=100
        ) # type: ignore 

    bbox_size: bpy.props.FloatProperty(
        name="BBOX Size", 
        default=1000,
        update=update_bbox_size  # Set the update callback
        )# type: ignore 

    # value_convert : bpy.props.EnumProperty(
    #     items=value_convert_items, 
    #     name="Value Type"
    #     )# type: ignore
    
    extracted_type : bpy.props.EnumProperty(
        items=extracted_type_items, 
        name="Extracted Type"
        )# type: ignore
    
    dense_type : bpy.props.EnumProperty(
        items=dense_type_items, 
        name="Dense Type"
        )# type: ignore
    
    dense_norm : bpy.props.EnumProperty(
        items=dense_norm_items, 
        name="Dense Norm"
        )# type: ignore

    particle_fix_size: bpy.props.FloatProperty(
        name="Particle Size", 
        min=0.0,
        max=100000.0,
        default=0.0
        )# type: ignore

    filter_min: bpy.props.FloatProperty(
        name="Min", 
        default=-1e+39
        )# type: ignore
    
    filter_max: bpy.props.FloatProperty(
        name="Max", 
        default=1e+39
        )# type: ignore
    
    density: bpy.props.FloatProperty(
        name="Density", 
        default=1.0
        )# type: ignore 

    anim_type : bpy.props.EnumProperty(
        items=anim_items, 
        name="Anim Type"
        )# type: ignore
    
    anim_frame: bpy.props.IntProperty(
        name="Frame", 
        default=0,
        min=0,
        max=1000000
        ) # type: ignore
    
    anim_start: bpy.props.IntProperty(
        name="Start Frame", 
        default=0,
        min=0,
        max=1000000
        ) # type: ignore

    anim_end: bpy.props.IntProperty(
        name="End Frame", 
        default=0,
        min=0,
        max=1000000
        ) # type: ignore
    
    anim_task_counter: bpy.props.IntProperty(
        name="Task Counter", 
        default=0,
        min=0,
        max=1000000
        ) # type: ignore    
    
    register_export: bpy.props.BoolProperty(
        name="Register Export", 
        default=False,
        update=lambda self, context: self.toggle_export_register(context)
        )# type: ignore

    slice_nr: bpy.props.IntProperty(
        name="Slice Number", 
        default=0,
        min=-1000000,
        max=1000000
        ) # type: ignore    

    slice_axis: bpy.props.EnumProperty(
        items=slice_axis_items,
        name="Slice Axis",        
        ) # type: ignore
    
    replace_path_enabled: bpy.props.BoolProperty(
        name="Replace Path",
        default=False,
        )# type: ignore
    
    replace_path_orig: bpy.props.StringProperty(
        name="Orig Path",
        default="",
        )# type: ignore
    
    replace_path_new: bpy.props.StringProperty(
        name="New Path",
        default="",
        )# type: ignore
    
    def enable_export_register(self, context):
        """Enable file export registration"""
        #if not self.register_export:
        #self.register_export = True
        print("Export registration enabled")
        self.anim_frame = context.scene.frame_current        
        # Add registration logic here
        if bpy.app.handlers.frame_change_post.count(self.frame_change_callback) == 0:
            bpy.app.handlers.frame_change_post.append(self.frame_change_callback)
    
    def disable_export_register(self, context):
        """Disable file export registration"""
        #if self.register_export:
        #    self.register_export = False
        print("Export registration disabled")
        # Add unregistration logic here
        if self.frame_change_callback in bpy.app.handlers.frame_change_post:
            bpy.app.handlers.frame_change_post.remove(self.frame_change_callback)
            
    def toggle_export_register(self, context):
        """Toggle export registration based on property value"""
        if self.register_export:
            # Registration logic
            self.enable_export_register(context)
        else:
            # Unregistration logic            
            self.disable_export_register(context)
        
    @staticmethod
    def frame_change_callback(scene, depsgraph):
        """Callback for frame change to extract data at current frame"""
        # Update frame in settings
        scene.view_pg_bspace.anim_frame = scene.frame_current
        # Extract data
        if scene.view_pg_bspace.register_export:
            scene.bspace.extract_data(bpy.context)
    

########################### OPENGL ###################################

def lerp_color(color1, color2, factor):
    """ Linearly interpolate between two colors """
    return [c1 * (1 - factor) + c2 * factor for c1, c2 in zip(color1, color2)]

def draw_callback_px():
    font_id = 0  # default font
    blf_size = 18.0
    blf.size(font_id, blf_size)  # set font size
    blf.color(font_id, 1.0, 1.0, 1.0, 1.0)  # Set font color to white

    x = 40
    y = 60
    width = 20
    height = 300  # total height for the legend
    num_steps = 100
    step_height = height / num_steps

    # start_color = (0.0, 0.0, 1.0, 1.0)  # Blue
    # mid_color = (0.0, 1.0, 0.0, 1.0)    # Green
    # end_color = (1.0, 0.0, 0.0, 1.0)    # Red

    obj_vdb = bpy.context.view_layer.objects.active
    if not obj_vdb is None and 'BSPACE' in obj_vdb:
        mat = obj_vdb.data.materials[0]
        
        #num_elements = len(mat.node_tree.nodes["Color Ramp"].color_ramp.elements)
        #mat.node_tree.nodes["Color Ramp"].color_ramp.elements[2].position
        #mat.node_tree.nodes["Color Ramp"].color_ramp.elements[2].color

        shader = gpu.shader.from_builtin('UNIFORM_COLOR')
        
        #step = int(num_steps / num_elements)

        # for j in range(num_elements - 1):
        #     c1 = mat.node_tree.nodes["Color Ramp"].color_ramp.elements[j + 0].color
        #     c2 = mat.node_tree.nodes["Color Ramp"].color_ramp.elements[j + 1].color
        #     p1 = mat.node_tree.nodes["Color Ramp"].color_ramp.elements[j + 0].position
        #     p2 = mat.node_tree.nodes["Color Ramp"].color_ramp.elements[j + 1].position

        #     step = int(num_steps / num_elements)
            #step = int(p1 * num_steps)

        for i in range(num_steps):
            #i = k + j * step
            # Interpolate between blue and green
            #factor = (k / step)
            #color = lerp_color(c1, c2, factor)            
            # # Interpolate between green and red
            # factor = ((i - num_steps / 2) / (num_steps / 2))
            # color = lerp_color(mid_color, end_color, factor)
            color = mat.node_tree.nodes["Color Ramp"].color_ramp.evaluate(i / num_steps)

            rect_y = y + i * step_height
            vertices = (
                (x, rect_y), (x + width, rect_y),
                (x + width, rect_y + step_height), (x, rect_y + step_height)
            )
            indices = ((0, 1, 2), (2, 3, 0))
            batch = batch_for_shader(shader, 'TRIS', {"pos": vertices}, indices=indices)
            shader.bind()
            shader.uniform_float("color", color)
            batch.draw(shader)    

        # Draw temperature labels
        blf.position(font_id, x + width + 5, y - blf_size / 2, 0)
        min_v = f'{obj_vdb["MIN_VALUE_REDUCED"]:.2E}'
        blf.draw(font_id, min_v)

        # blf.position(font_id, x + width + 5, y + height / 2 - 10, 0)
        # blf.draw(font_id, "10e5")
        blf.position(font_id, x + width + 5, y + height - blf_size / 2, 0)
        max_v = f'{obj_vdb["MAX_VALUE_REDUCED"]:.2E}'
        blf.draw(font_id, max_v)
        
        # Draw labels
        blf.position(font_id, x, y + height + blf_size / 2, 0)
        #blf.draw(font_id, "GASS - RHO")    
        desc = obj_vdb['NAME'] #str(particle_type_items[int(obj_vdb['PARTICLE_TYPE'])][1]) + str(" - ") + str(block_name_items[int(obj_vdb['BLOCK_NAME'])][1])
        blf.draw(font_id, desc)

########################### utils #############################
def get_view_bounds_3d(context):
    # This function assumes we're dealing with the 3D view

    # Get the current region and region data
    #region = context.region
    region_3d = context.space_data.region_3d

    if region_3d.view_perspective == 'ORTHO':
        # Get the perspective matrix
        perspective_matrix = region_3d.perspective_matrix

        # Inverse perspective matrix to transform from clip space to world space
        inv_perspective_matrix = perspective_matrix.inverted()

        # Define clip space coordinates
        clip_coords = [mathutils.Vector((-1.0, -1.0, -1.0)),
                    mathutils.Vector((1.0, -1.0, -1.0)),
                    mathutils.Vector((1.0, 1.0, -1.0)),
                    mathutils.Vector((-1.0, 1.0, -1.0)),
                    mathutils.Vector((-1.0, -1.0, 1.0)),
                    mathutils.Vector((1.0, -1.0, 1.0)),
                    mathutils.Vector((1.0, 1.0, 1.0)),
                    mathutils.Vector((-1.0, 1.0, 1.0))]

        # Transform clip space coordinates to world space
        world_coords = [inv_perspective_matrix @ coord for coord in clip_coords]

        # Use numpy for vectorized min/max calculation
        coords_np = np.array([[v.x, v.y, v.z] for v in world_coords])

        # Calculate min/max for the bounding box in world space
        min_bounds = coords_np.min(axis=0)
        max_bounds = coords_np.max(axis=0)

        #min_abs_value = min(tuple(min_bounds) + tuple(max_bounds), key=abs)

        grid_dim = context.scene.view_pg_bspace.grid_dim * float(context.scene.view_pg_bspace.bbox_size) / float(region_3d.view_distance)

        if grid_dim < context.scene.view_pg_bspace.grid_dim:
            grid_dim = context.scene.view_pg_bspace.grid_dim
    
    else:
        min_bounds = (0,0,0)
        max_bounds = (context.scene.view_pg_bspace.bbox_size,context.scene.view_pg_bspace.bbox_size,context.scene.view_pg_bspace.bbox_size)

        grid_dim = context.scene.view_pg_bspace.grid_dim

    #print(f"Bounding Box Min: {mathutils.Vector(min_bounds)}, Max: {mathutils.Vector(max_bounds)}, GridDim: {int(grid_dim)}")

    return mathutils.Vector(min_bounds), mathutils.Vector(max_bounds) #, int(grid_dim)

def calculate_view_bounding_box_from_obj(context):
    # Ensure we are in object mode
    bpy.ops.object.mode_set(mode='OBJECT')

    # Get the current view matrix
    view_matrix = context.space_data.region_3d.view_matrix
    # Invert the view matrix to transform vertices to view space
    inverse_view_matrix = view_matrix.inverted()

    # Initialize min and max vectors
    min_vec = mathutils.Vector((float('inf'), float('inf'), float('inf')))
    max_vec = mathutils.Vector((float('-inf'), float('-inf'), float('-inf')))

    # Loop through all visible objects in the scene
    for obj in context.visible_objects:
        # Make sure the object has geometry data
        if obj.type == 'MESH':
            # Apply object's matrix_world to get vertices in world space, then apply inverse_view_matrix
            world_vertices = [inverse_view_matrix @ obj.matrix_world @ v.co for v in obj.data.vertices]

            # Update min and max vectors
            for vert in world_vertices:
                min_vec = mathutils.Vector((min(min_vec.x, vert.x), min(min_vec.y, vert.y), min(min_vec.z, vert.z)))
                max_vec = mathutils.Vector((max(max_vec.x, vert.x), max(max_vec.y, vert.y), max(max_vec.z, vert.z)))

    # Now min_vec and max_vec define the bounding box in view space
    #print(f"Bounding Box Min: {min_vec}, Max: {max_vec}")

    return min_vec, max_vec

def toggle_perspective_ortho(context):
    for area in context.screen.areas:
        if area.type == 'VIEW_3D':
            # Get the space data from the 3D view area
            space_data = area.spaces.active
            # Check the current viewport perspective and toggle
            if space_data.region_3d.view_perspective == 'PERSP':
                space_data.region_3d.view_perspective = 'ORTHO'
            else:
                space_data.region_3d.view_perspective = 'PERSP'
            break


######################## OPERATORS ######################################
class BSPACE_OT_move_to_click(bpy.types.Operator):
    """Move Object to Clicked Point"""
    bl_idname = "bspace.move_to_click"
    bl_label = "Move BBOX to Click"
    bl_options = {'REGISTER', 'UNDO'}

    def modal(self, context, event):
        if event.type == 'MOUSEMOVE':
            self.mouse_path = (event.mouse_region_x, event.mouse_region_y)
        
        elif event.type == 'LEFTMOUSE':
            # Get the mouse coordinates
            mouse_coords = (event.mouse_region_x, event.mouse_region_y)
            # Convert mouse coordinates to 3D space coordinates
            location = view3d_utils.region_2d_to_location_3d(
                context.region, context.region_data, mouse_coords, (0, 0, 0))
            
            if location:
                # Move the active object to this location
                # context.active_object.location = location
                context.scene.bspace.select_bbox(context)
                bbox = context.scene.bspace.get_bbox(context)
                bbox.location = location

                context.area.tag_redraw()
                return {'FINISHED'}
        
        elif event.type in {'RIGHTMOUSE', 'ESC'}:
            return {'CANCELLED'}

        return {'RUNNING_MODAL'}

    def invoke(self, context, event):
        if context.object:
            self.mouse_path = (0,0)
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "No active object, could not finish")
            return {'CANCELLED'}
        
class BSPACE_OT_move_to_cursor(bpy.types.Operator):
    bl_idname = "bspace.move_to_cursor"
    bl_label = "Move BBOX to Cursor"

    def execute(self, context):
        context.scene.bspace.select_bbox(context)
        bbox = context.scene.bspace.get_bbox(context)        
        bbox.location = context.scene.cursor.location

        return {'FINISHED'}        
        
class BSPACE_OT_LoadData(bpy.types.Operator):
    bl_idname = "bspace.load_data"
    bl_label = "Load Data"

    def execute(self, context):
        context.scene.bspace.load_data(context)
        return {'FINISHED'}
    
# class BSPACE_OT_Zoom(bpy.types.Operator):
#     bl_idname = "bspace.zoom"
#     bl_label = "Zoom"

#     def execute(self, context):
#         #context.scene.bspace.connect(context)
#         #context.scene.bspace.find_particle_data_types(context)
#         bpy.ops.view3d.zoom_border('INVOKE_DEFAULT')

#         print("bspace.zoom")
#         return {'FINISHED'}

class BSPACE_OT_Zoom(bpy.types.Operator):
    """Operator that zooms to the border and notifies on completion"""
    bl_idname = "bspace.zoom"
    bl_label = "Zoom Border and Notify"

    _initial_view_matrix = None
    
    def modal(self, context, event):
        # if event.type in {'LEFTMOUSE', 'RIGHTMOUSE', 'ESC'}:
        #     # Print message when zoom operation is actually completed or cancelled
        #     print("Zoom border operation completed or cancelled.")
        #     return {'FINISHED'}
        
        current_view_matrix = context.space_data.region_3d.view_matrix
        if current_view_matrix != self._initial_view_matrix:
            print("Viewport change detected after zoom operation.")
            context.scene.bspace.extract_data(context)
            return {'FINISHED'}
            #return {'PASS_THROUGH'}
                
        return {'PASS_THROUGH'}

    def invoke(self, context, event):
        # Store the initial view matrix before starting the zoom operation
        self._initial_view_matrix = context.space_data.region_3d.view_matrix.copy()

        # Initialize zoom border operation
        context.window_manager.modal_handler_add(self)
        bpy.ops.view3d.zoom_border('INVOKE_DEFAULT')
        return {'RUNNING_MODAL'}

class BSPACE_OT_Connect(bpy.types.Operator):
    bl_idname = "bspace.connect"
    bl_label = "Connect"

    def execute(self, context):
        context.scene.bspace.connect(context)
        context.scene.bspace.find_particle_data_types(context)        
        return {'FINISHED'}


class BSPACE_OT_Disconnect(bpy.types.Operator):
    bl_idname = "bspace.disconnect"
    bl_label = "Disconnect"

    def execute(self, context):
        context.scene.bspace.disconnect(context)
        return {'FINISHED'}
    
class BSPACE_OT_ExtractData(bpy.types.Operator):
    bl_idname = "bspace.extract_data"
    bl_label = "Extract Data"

    def execute(self, context):
        try:
            context.scene.bspace.extract_data(context)
        except:
            #reconnect
            context.scene.bspace.disconnect(context)
            context.scene.bspace.connect(context)

            #try again
            context.scene.bspace.extract_data(context)

        return {'FINISHED'}
    
class BSPACE_OT_SelectBBOX(bpy.types.Operator):
    bl_idname = "bspace.select_bbox"
    bl_label = "Select BBOX"

    def execute(self, context):
        context.scene.bspace.select_bbox(context)
        return {'FINISHED'}

class BSPACE_OT_FindBBOX(bpy.types.Operator):
    bl_idname = "bspace.find_bbox"
    bl_label = "Find BBOX"

    def execute(self, context):
        context.scene.bspace.find_bbox(context)
        return {'FINISHED'}     

class BSPACE_OT_CreateBBox(bpy.types.Operator):
    bl_idname = "bspace.create_bbox"
    bl_label = "Create BBOX"

    def execute(self, context):
        context.scene.bspace.create_bbox(context)
        return {'FINISHED'}
    
class BSPACE_OT_FindTypes(bpy.types.Operator):
    bl_idname = "bspace.find_types"
    bl_label = "Find Types"

    def execute(self, context):
        context.scene.bspace.find_particle_data_types(context)
        return {'FINISHED'}    


class BSPACE_OT_SetWorldEnvironment(bpy.types.Operator):
    bl_idname = "bspace.set_world_environment"
    bl_label = "Set World Environment"

    def execute(self, context):
        context.scene.bspace.set_world_environment(context)
        context.scene.bspace.set_view(context)
        return {'FINISHED'}     
    
class BSPACE_OT_select(bpy.types.Operator):
    """Select the object from the list."""
    bl_idname = "bspace.select_from_list"
    bl_label = "Select Object"

    index: bpy.props.IntProperty() # type: ignore

    def execute(self, context):
        
        # hide all vdbs
        for ob in bpy.data.objects:
            if 'BSPACE' in ob:
                ob.hide_render = True
                ob.hide_set(True)

        ob = bpy.data.objects[context.scene.bspace_list_data[self.index].Obj]
        context.view_layer.objects.active = ob
        ob.hide_render = False
        ob.hide_set(False)

        return {'FINISHED'}
    
class BSPACE_OT_CalcBBOXViewport(bpy.types.Operator):
    bl_idname = "bspace.calc_bbox_viewport"
    bl_label = "Calc BBOX Viewport"

    def execute(self, context):
        get_view_bounds_3d(context)
        #bpy.data.objects["Cube"].location = ((max_vec[0]-min_vec[0])/2.0,(max_vec[1]-min_vec[1])/2.0,(max_vec[2]-min_vec[2])/2.0)

        return {'FINISHED'}
    
class BSPACE_OT_MergeVDB(bpy.types.Operator):
    bl_idname = "bspace.merge_vdb"
    bl_label = "MergeVDB"

    def execute(self, context):
        context.scene.bspace.merge_vdb(context)

        return {'FINISHED'}
    
class BSPACE_OT_HistoVDB(bpy.types.Operator):
    bl_idname = "bspace.vdb2histo"
    bl_label = "Histogram VDB"

    def execute(self, context):
        context.scene.bspace.vdb2histo(context)

        return {'FINISHED'}

class BSPACE_OT_VDB2PNG(bpy.types.Operator):
    bl_idname = "bspace.vdb2png"
    bl_label = "VDB to PNG"

    def execute(self, context):
        context.scene.bspace.vdb2png(context)

        return {'FINISHED'}

class BSPACE_OT_MergeDeselectAll(bpy.types.Operator):
    bl_idname = "bspace.merge_deselect_all"
    bl_label = "Deselect All"

    def execute(self, context):
        # Deselect all objects first
        bpy.ops.object.select_all(action='DESELECT')

        #pref = haystack_pref.preferences()
        # hide all vdbs
        for ob in bpy.data.objects:
            if 'BSPACE' in ob:
                ob.hide_render = True
                ob.hide_set(True)

        return {'FINISHED'}     
##############################################################

class BSPACE_OT_update_extracted_data(bpy.types.Operator):
    bl_idname = 'bspace.update_extracted_data'
    bl_label = 'Update'

    name : bpy.props.StringProperty(
        ) # type: ignore
    
    id : bpy.props.IntProperty(
        ) # type: ignore

    def execute(self, context):

        # Deselect all objects first
        bpy.ops.object.select_all(action='DESELECT')

        #pref = haystack_pref.preferences()
        # hide all vdbs
        for ob in bpy.data.objects:
            if 'BSPACE' in ob:
                ob.hide_render = True
                ob.hide_set(True)

        ob = bpy.data.objects[self.name]
        context.view_layer.objects.active = ob
        ob.hide_render = False
        ob.hide_set(False)
        ob.select_set(True)

        return {'FINISHED'}

class BSPACE_OT_merge_select_data(bpy.types.Operator):
    bl_idname = 'bspace.merge_select_data'
    bl_label = 'Select'

    name : bpy.props.StringProperty(
        ) # type: ignore
    
    id : bpy.props.IntProperty(
        ) # type: ignore

    def execute(self, context):

        # # Deselect all objects first
        # bpy.ops.object.select_all(action='DESELECT')

        # #pref = haystack_pref.preferences()
        # # hide all vdbs
        # for ob in bpy.data.objects:
        #     if 'BSPACE' in ob:
        #         ob.hide_render = True
        #         ob.hide_set(True)

        ob = bpy.data.objects[self.name]
        context.view_layer.objects.active = ob
        ob.hide_render = False
        ob.hide_set(False)
        ob.select_set(True)

        return {'FINISHED'}         
    
class BSPACE_PG_ExtractedData(bpy.types.PropertyGroup):
    Id : bpy.props.IntProperty(name="Id") # type: ignore
    Name : bpy.props.StringProperty(name="Name") # type: ignore
    Obj : bpy.props.StringProperty() # type: ignore

class BSPACE_UL_ExtractedData(bpy.types.UIList):
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname):
        op = layout.operator("bspace.update_extracted_data", text=item.Name, icon='FILE_VOLUME')
        op.name = item.Name
        op.id = item.Id

    #     if item:
    #         layout.label(text=('%d' % item.Id))
    #         layout.label(text=item.Name)

    # def filter_items(self, context, data, propname):
    #     """Filter and order items in the list."""

    #     filtered = []
    #     ordered = []

    #     items = getattr(data, propname)

    #     helpers = bpy.types.UI_UL_list
    #     filtered = helpers.filter_items_by_name(self.filter_name,
    #                                     self.bitflag_filter_item,
    #                                     items, "Name", reverse=False)

    #     return filtered, ordered

class BSPACE_UL_MergeData(bpy.types.UIList):
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname):
        op = layout.operator("bspace.merge_select_data", text=item.Name, icon='FILE_VOLUME')
        op.name = item.Name
        op.id = item.Id

class BSPACE_PG_ExtractedTypes(bpy.types.PropertyGroup):
    Id : bpy.props.IntProperty(name="Id") # type: ignore
    Name : bpy.props.StringProperty(name="Name") # type: ignore

    Type : bpy.props.StringProperty(name="Type") # type: ignore
    Type_id : bpy.props.IntProperty(name="Type_id") # type: ignore
    Block : bpy.props.StringProperty(name="Block") # type: ignore
    Block_id : bpy.props.IntProperty(name="Block_id") # type: ignore

class BSPACE_UL_ExtractedTypes(bpy.types.UIList):
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname):
        if item:
            #layout.label(text=('%d' % item.Id))
            layout.label(text=item.Name)

    def filter_items(self, context, data, propname):
        """Filter and order items in the list."""

        filtered = []
        ordered = []

        items = getattr(data, propname)

        helpers = bpy.types.UI_UL_list
        filtered = helpers.filter_items_by_name(self.filter_name,
                                        self.bitflag_filter_item,
                                        items, "Name", reverse=False)

        return filtered, ordered                   

##############################################################
class BSPACE:
    def __init__(self):
        self.enabled = False

        #self.connection_thread = None
        #self.connection_ssh1 = None
        #self.connection_ssh2 = None
        #self.connection_forward_server = None
        self.connection_channel = None

        self.bspace_process = None
        self.bspace_tunnel = None

        self.draw_handler = None

    def set_world_environment(self, context):
        bpy.data.worlds["World"].use_nodes = True
        bpy.data.worlds["World"].node_tree.nodes["Background"].inputs[0].default_value = (0,0,0,1)
        bpy.data.worlds["World"].node_tree.nodes["Background"].inputs[1].default_value = 0

        bpy.data.worlds["World"].color = (0,0,0)

        shading = context.area.spaces.active.shading
        shading.background_type = 'WORLD'
        shading.type = 'RENDERED' #'MATERIAL'
        shading.use_scene_world = True

        context.scene.eevee.volumetric_start = 0.001
        context.scene.eevee.volumetric_end = 1000.0
        context.scene.eevee.volumetric_tile_size = '4'
        context.scene.eevee.volumetric_samples = 64

        context.scene.render.engine = 'CYCLES'
        context.scene.cycles.device = 'GPU'

        context.scene.cycles.volume_preview_step_rate = 0.1
        context.scene.cycles.volume_max_steps = 64

        pass

    def set_vdb_shader(self, context, obj_vdb):
        if obj_vdb is None:
            obj_vdb = context.view_layer.objects.active

        obj_vdb.data.materials.clear()

        # Create a new material
        mat = bpy.data.materials.new(name="VolumeMat")
        mat.use_nodes = True
        nodes = mat.node_tree.nodes

        # Clear default nodes
        for node in nodes:
            nodes.remove(node)

        # Create an Attribute node
        attribute_node = nodes.new(type='ShaderNodeAttribute')
        attribute_node.attribute_name = "density"  # Set the attribute name to 'density'
        attribute_node.location = (-200, 0)

        maprange_node = nodes.new(type='ShaderNodeMapRange')
        #attribute_node.attribute_name = "density"  # Set the attribute name to 'density'
        # obj_vdb_new['MIN_VALUE'] = min_value
        # obj_vdb_new['MAX_VALUE'] = max_value
        maprange_node.inputs['From Min'].default_value = obj_vdb['MIN_VALUE_REDUCED']
        maprange_node.inputs['From Max'].default_value = obj_vdb['MAX_VALUE_REDUCED']
        maprange_node.location = (0, 0)

        # Create a ColorRamp node
        color_ramp_node = nodes.new(type='ShaderNodeValToRGB')
        color_ramp_node.color_ramp.elements[0].color = (0, 0, 1, 0)  # Blue at position 0
        color_ramp_node.color_ramp.elements[1].color = (1, 0, 0, 1)  # Red at position 1
        color_ramp_node.color_ramp.elements.new(0.5)  # Add a new stop at 0.5
        color_ramp_node.color_ramp.elements[1].color = (0, 1, 0, 0.5)  # Change this new stop to Green
        color_ramp_node.location = (200, 0)

        # Create a Math node
        math_node = nodes.new(type='ShaderNodeMath')
        math_node.operation = 'MULTIPLY'
        math_node.inputs[0].default_value = 1.0
        math_node.inputs[1].default_value = 1.0
        math_node.location = (500, -60)

        # Create a Emission node
        emission_node = nodes.new(type='ShaderNodeEmission')
        emission_node.location = (700, 0)

        # Create a VolumeAbsorption node
        volume_absorption_node = nodes.new(type='ShaderNodeVolumeAbsorption')
        volume_absorption_node.location = (700, -200)   
        volume_absorption_node.inputs['Color'].default_value = (0.3, 0.3, 0.3, 1)

        # Create a AddShader node
        add_shader_node = nodes.new(type='ShaderNodeAddShader')
        add_shader_node.location = (1000, 0)   

        # Create an Output node
        output_node = nodes.new(type='ShaderNodeOutputMaterial')
        output_node.location = (1300, 0)

        # Link nodes
        links = mat.node_tree.links
        links.new(attribute_node.outputs['Fac'], maprange_node.inputs['Value'])
        links.new(maprange_node.outputs['Result'], color_ramp_node.inputs['Fac'])        
        
        links.new(color_ramp_node.outputs['Color'], emission_node.inputs['Color'])

        links.new(color_ramp_node.outputs['Alpha'], math_node.inputs[0])
        links.new(math_node.outputs['Value'], emission_node.inputs['Strength'])        

        links.new(color_ramp_node.outputs['Alpha'], volume_absorption_node.inputs['Density'])

        links.new(emission_node.outputs['Emission'], add_shader_node.inputs[0])
        links.new(volume_absorption_node.outputs['Volume'], add_shader_node.inputs[1])

        links.new(add_shader_node.outputs['Shader'], output_node.inputs['Volume'])

        # Assign the material to the active object
        obj_vdb.data.materials.append(mat)

        try:
            context.scene.haystack.server_settings.mat_volume = mat
        except:
            pass

        #############################################

        # Remember to remove the draw handler when it's no longer needed
        if self.draw_handler:
            bpy.types.SpaceView3D.draw_handler_remove(self.draw_handler, 'WINDOW')

        # Add draw handler
        self.draw_handler = bpy.types.SpaceView3D.draw_handler_add(draw_callback_px, (), 'WINDOW', 'POST_PIXEL')

        pass

    def set_vdb_shader2(self, context, obj_vdb):
        if obj_vdb is None:
            obj_vdb = context.view_layer.objects.active

        obj_vdb.data.materials.clear()

        # Create a new material
        mat = bpy.data.materials.new(name="VolumeMaterial")
        mat.use_nodes = True
        nodes = mat.node_tree.nodes

        # Clear default nodes
        for node in nodes:
            nodes.remove(node)

        # Create an Attribute node
        attribute_node = nodes.new(type='ShaderNodeAttribute')
        attribute_node.attribute_name = "density"  # Set the attribute name to 'density'
        attribute_node.location = (0, 0)

        # Create a ColorRamp node
        color_ramp_node = nodes.new(type='ShaderNodeValToRGB')
        color_ramp_node.color_ramp.elements[0].color = (0, 0, 1, 0)  # Blue at position 0
        color_ramp_node.color_ramp.elements[1].color = (1, 0, 0, 1)  # Red at position 1
        color_ramp_node.color_ramp.elements.new(0.5)  # Add a new stop at 0.5
        color_ramp_node.color_ramp.elements[1].color = (0, 1, 0, 0.5)  # Change this new stop to Green
        color_ramp_node.location = (200, 0)

        # Create a Principled Volume node
        principled_volume_node = nodes.new(type='ShaderNodeVolumePrincipled')
        principled_volume_node.inputs['Density'].default_value = 1.0  # Default density
        principled_volume_node.inputs['Emission Strength'].default_value = 1.0  # Default emission strength
        principled_volume_node.inputs['Density Attribute'].default_value = ''
        principled_volume_node.location = (600, 0)

        # Create an Output node
        output_node = nodes.new(type='ShaderNodeOutputMaterial')
        output_node.location = (1000, 0)

        # Link nodes
        links = mat.node_tree.links
        links.new(attribute_node.outputs['Fac'], color_ramp_node.inputs['Fac'])
        
        links.new(color_ramp_node.outputs['Color'], principled_volume_node.inputs['Color'])
        links.new(color_ramp_node.outputs['Color'], principled_volume_node.inputs['Emission Color'])
        links.new(color_ramp_node.outputs['Alpha'], principled_volume_node.inputs['Density'])
        links.new(color_ramp_node.outputs['Alpha'], principled_volume_node.inputs['Emission Strength'])

        links.new(principled_volume_node.outputs['Volume'], output_node.inputs['Volume'])

        # Assign the material to the active object
        obj_vdb.data.materials.append(mat)


        #############################################

        # Remember to remove the draw handler when it's no longer needed
        if self.draw_handler:
            bpy.types.SpaceView3D.draw_handler_remove(self.draw_handler, 'WINDOW')

        # Add draw handler
        self.draw_handler = bpy.types.SpaceView3D.draw_handler_add(draw_callback_px, (), 'WINDOW', 'POST_PIXEL')

        pass   

    def set_view(self, context):
        space_data = context.area.spaces.active
        
        space_data.clip_start = 0.1
        space_data.clip_end = 100000

        pass         

    def create_instances_bobj(self, context, obj_collection, name):
        empty_name = '%s_instances' % name
        mesh = bpy.data.meshes.new(name=empty_name)

        mesh.update()
        mesh.validate()

        obj = bpy.data.objects.new(empty_name, mesh)
        obj_collection.objects.link(obj)

        return obj

    def create_instances_gnode(self, context, obj, instance_object):
        # pass
        # context.view_layer.objects.active = obj
        bpy.ops.node.new_geometry_nodes_modifier()

        obj.modifiers[0].name = instance_object.name
        node_group = obj.modifiers[0].node_group
        node_group.name = instance_object.name

        # bpy.ops.node.add_node(type="GeometryNodeObjectInfo", use_transform=True)
        node_object_info = node_group.nodes.new('GeometryNodeObjectInfo')
        node_object_info.location = (-400, 600)
        node_object_info.inputs[0].default_value = instance_object
        node_object_info.inputs[1].default_value = True
        # bpy.ops.node.add_node(type="GeometryNodeInstanceOnPoints", use_transform=True)
        node_instance = node_group.nodes.new('GeometryNodeInstanceOnPoints')
        node_instance.location = (-50, 100)

        # bpy.ops.node.add_node(type="GeometryNodeInputNamedAttribute", use_transform=True)
        # node_in_rot = node_group.nodes.new('GeometryNodeInputNamedAttribute')
        # node_in_rot.data_type = 'FLOAT_VECTOR'
        # node_in_rot.inputs[0].default_value = 'rot'
        # node_in_rot.location = (-400, 350)
        # bpy.ops.node.add_node(type="GeometryNodeInputNamedAttribute", use_transform=True)
        node_in_scal = node_group.nodes.new('GeometryNodeInputNamedAttribute')
        node_in_scal.data_type = 'FLOAT'
        node_in_scal.inputs[0].default_value = 'value'
        node_in_scal.location = (-400, 200)

        node_in = node_group.nodes['Group Input']
        node_out = node_group.nodes['Group Output']

        node_group.links.clear()

        node_group.links.new(
            node_in.outputs['Geometry'], node_instance.inputs['Points'])
        # node_group.links.new(
        #     node_in_rot.outputs['Attribute'], node_instance.inputs['Rotation'])
        # node_group.links.new(
        #     node_in_scal.outputs['Attribute'], node_instance.inputs['Scale'])

        node_group.links.new(
            node_object_info.outputs['Geometry'], node_instance.inputs['Instance'])

        node_group.links.new(
            node_instance.outputs['Instances'], node_out.inputs['Geometry'])

    def add_verts_to_mesh(self, context, empty_obj, count):
        context.view_layer.objects.active = empty_obj

        mesh = empty_obj.data

        mesh.vertices.add(count)

        # mesh.attributes.new(name="loc", type='FLOAT_VECTOR', domain='POINT')
        #mesh.attributes.new(name="rot", type='FLOAT_VECTOR', domain='POINT')
        mesh.attributes.new(name="value", type='FLOAT', domain='POINT')

        mesh.update()
        mesh.validate()                    

    def load_data(self, context):
        # for x in os.listdir(context.scene.view_pg_bspace.local_data_dir_path):
        #     if x.endswith(".vdb"):
        #         # Prints only text file present in My Folder
        #         print(x)

        #     if x.endswith(".bin"):
        #         # Prints only text file present in My Folder
        #         print(x)
        object_collection = bpy.data.collections.new('bspace')
        context.scene.collection.children.link(object_collection)

        for x in pathlib.Path(context.scene.view_pg_bspace.local_data_dir_path).glob("*.vdb"):    
            print(x)             
            bpy.ops.object.volume_import(filepath=str(x), align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
            volume_obj = bpy.context.active_object
            volume_obj.hide_render = True
            volume_obj.hide_set(True)            
        
        import struct
        for x in pathlib.Path(context.scene.view_pg_bspace.local_data_dir_path).glob("*.bin"):    
            print(x)                 

            obj = bpy.data.objects['instance']
            with open(str(x), "rb") as f:
                # Read no_points of type size_t. Assuming it's an 8-byte unsigned integer.
                # This assumption is based on typical 64-bit architectures.
                no_values = struct.unpack('Q', f.read(8))[0]
                no_pos = no_values * 3

                # Read positions
                positions = [struct.unpack('f', f.read(4))[0] for _ in range(no_pos)]

                # Read values
                values = [struct.unpack('f', f.read(4))[0] for _ in range(no_values)]

                #return no_points, positions, values

                #return dimensions, spacing, pixels                                            
                empty_obj = self.create_instances_bobj(context, object_collection, x.name)
                empty_obj.hide_render = True
                empty_obj.hide_set(True)

                self.add_verts_to_mesh(context, empty_obj, no_values)

                #for i in range(no_values):
                #    context.object.data.vertices[i].co = (positions[0 + i * 3], positions[1 + i * 3], positions[2 + i * 3])
                #    context.object.data.attributes["value"].data[i].value = values[i]
                context.object.data.vertices.foreach_set("co", positions)
                context.object.data.attributes["value"].data.foreach_set("value", values)

                self.create_instances_gnode(context, empty_obj, obj)        

        pass

    def int_to_bytes(self, i):
        return struct.pack('i', int(i))
    
    def int3_to_bytes(self, i):
        return struct.pack('iii', int(i[0]), int(i[1]), int(i[2]))
    
    def float_to_bytes(self, f):
        return struct.pack('f', float(f))
    
    def float3_to_bytes(self, f):
        return struct.pack('fff', float(f[0]), float(f[1]), float(f[2]))
    
    def bytes_to_int(self, i):
        return struct.unpack('i', i)[0]
    
    def bytes_to_int64(self, i):
        return struct.unpack('q', i)[0]
    
    def bytes_to_float(self, f):
        return struct.unpack('f', f)[0]
    
    def bytes_to_double(self, d):
        return struct.unpack('d', d)[0]
    
    def bytes_to_float3(self, f):
        return struct.unpack('fff', f)
    
    def find_particle_data_types(self, context):
        #send
        message_type = 1
        bmessage_type = self.int_to_bytes(message_type)
        self.tcp_send(context, bmessage_type)

        # recv
        context.scene.view_pg_bspace.anim_type = str(self.bytes_to_int(self.tcp_recv(context, 4)))
        context.scene.view_pg_bspace.anim_start = self.bytes_to_int(self.tcp_recv(context, 4))
        context.scene.view_pg_bspace.anim_end = self.bytes_to_int(self.tcp_recv(context, 4))
        
		#tcp::send_data((char*)&s, sizeof(int));
        data_size = self.bytes_to_int(self.tcp_recv(context, 4))
		#tcp::send_data((char*)particle_data_types.c_str(), sizeof(char) * s);
        data_bytes = self.tcp_recv(context, data_size)        
        
        # tcp::recv_data((char*)&ack, sizeof(ack));
        self.tcp_send(context, bmessage_type)

        # Decode bytes to string
        data_str = data_bytes.decode('utf-8')

        # Initialize an array to hold the structures
        #structures = []

        # Split the string into lines
        entries = data_str.strip().split('\n')

        context.scene.bspace_list_types.clear()
        #context.scene.bspace_list_types_index = -1

        # Parse each line
        for entry in entries:
            fields = entry.split(';')
            if len(fields) == 4:
                # Map fields to a dictionary (or a similar structure if needed)
                # structure = {
                #     'type': fields[0],
                #     'type_id': int(fields[1]),
                #     'block': fields[2],
                #     'block_id': int(fields[3])
                # }
                #structures.append(structure)

                # add to list
                item = context.scene.bspace_list_types.add()
                item.Id = len(context.scene.bspace_list_types)                

                item.Type = fields[0]
                item.Type_id = int(fields[1])
                item.Block = fields[2]
                item.Block_id = int(fields[3])

                item.Name = item.Type + " - " + item.Block

        if context.scene.bspace_list_types_index > len(context.scene.bspace_list_types) - 1:
            context.scene.bspace_list_types_index = len(context.scene.bspace_list_types) - 1                 

        #return structures

    def get_bbox_name(self):
        return "BSPACE_BBOX"
    
    def set_bbox_size(self, context, value):
        bbox = self.get_bbox(context)

        #bbox.empty_display_size = value / 2
        bbox.scale=(value / 2, value / 2, value / 2)

    # def get_bbox_size(self, context):
    #     bbox = self.get_bbox(context)
    #     #return bbox.empty_display_size * 2
    #     return max(bbox.scale) * 2
    
    def get_bbox_dim(self, context):
        bbox = self.get_bbox(context)

        #return bbox.scale * bbox.empty_display_size * 2
        return bbox.scale * 2

    def get_bbox_location(self, context):
        bbox = self.get_bbox(context)

        return bbox.location

    def get_bbox(self, context):
        if not self.get_bbox_name() in bpy.data.objects:
            self.create_bbox(context)

        return bpy.data.objects[self.get_bbox_name()]

    def create_bbox(self, context):
        if self.get_bbox_name() in bpy.data.objects:
            obj = bpy.data.objects[self.get_bbox_name()]

            # If the object is linked to a collection, unlink it first
            for collection in obj.users_collection:
                collection.objects.unlink(obj)

            # Delete the object
            bpy.data.objects.remove(obj)

        bbox_size = context.scene.view_pg_bspace.bbox_size
        bbox_size_half = bbox_size / 2.0

        bpy.ops.object.empty_add(type='CUBE')        
        obj_select = bpy.context.view_layer.objects.active
        #context.scene.view_pg_bspace.box_object_select = obj_select

        #obj_select.empty_display_size = bbox_size_half
        #obj_select.empty_display_size = 1
        obj_select.scale=(bbox_size_half, bbox_size_half, bbox_size_half)
        obj_select.location=(bbox_size_half, bbox_size_half, bbox_size_half)

        obj_select.name = self.get_bbox_name()   

    def select_bbox(self, context):
        bbox = self.get_bbox(context)

        # Deselect all objects first
        bpy.ops.object.select_all(action='DESELECT')
        context.view_layer.objects.active = bbox        
        bbox.select_set(True)

    def find_bbox(self, context):
        # bbox = self.get_bbox(context)

        # # Deselect all objects first
        # bpy.ops.object.select_all(action='DESELECT')
        # context.view_layer.objects.active = bbox        
        # bbox.select_set(True)

        #     #send
        #     message_type = 1
        #     bmessage_type = self.int_to_bytes(message_type)
        #     self.tcp_send(context, bmessage_type)

        #     # tcp::recv_data((char*)&particle_type, sizeof(int));
        #     bparticle_type = self.int_to_bytes(int(context.scene.view_pg_bspace.particle_type))        
        #     self.tcp_send(context, bparticle_type)          

        #     # tcp::recv_data((char*)&block_name_id, sizeof(int));
        #     bblock_name_id = self.int_to_bytes(int(context.scene.view_pg_bspace.block_name))        
        #     self.tcp_send(context, bblock_name_id)        

        #     #recv
        #     snapshot_size = float(self.bytes_to_double(self.tcp_recv(context, 8))) # TODO
        #     #context.scene.view_pg_bspace.snapshot_size

        if context.scene.bspace_list_types_index < 0 or context.scene.bspace_list_types_index > len(context.scene.bspace_list_types):
            raise Exception("context.scene.bspace_list_types_index < 0 or context.scene.bspace_list_types_index > len(context.scene.bspace_list_types)")

        # bbox_location = self.get_bbox_location(context)
        # bbox_dims = self.get_bbox_dim(context)
        # bbox_size = self.get_bbox_size(context)

        # bbox_min = (bbox_location[0] - bbox_dims[0] / 2.0, bbox_location[1] - bbox_dims[1] / 2.0, bbox_location[2] - bbox_dims[2] / 2.0)
        # bbox_max = (bbox_location[0] + bbox_dims[0] / 2.0, bbox_location[1] + bbox_dims[1] / 2.0, bbox_location[2] + bbox_dims[2] / 2.0)       
        
        # grid_dim = context.scene.view_pg_bspace.grid_dim               

        message_type = 3 #eBBOX
        bmessage_type = self.int_to_bytes(message_type)
        self.tcp_send(context, bmessage_type)        

        #################send
        # tcp::recv_data((char*)&particle_type, sizeof(int));
        particle_type = context.scene.bspace_list_types[context.scene.bspace_list_types_index].Type_id
        bparticle_type = self.int_to_bytes(int(particle_type))        
        self.tcp_send(context, bparticle_type)          

        # tcp::recv_data((char*)&block_name_id, sizeof(int));
        block_name = context.scene.bspace_list_types[context.scene.bspace_list_types_index].Block_id
        bblock_name_id = self.int_to_bytes(int(block_name))        
        self.tcp_send(context, bblock_name_id)
        #################recv
        #tcp::recv_data((char*)&bbox_min[0], sizeof(float) * 3);
        bbbox_min = self.bytes_to_float3(self.tcp_recv(context, 4 * 3))     

        #tcp::recv_data((char*)&bbox_max[0], sizeof(float) * 3);
        bbbox_max = self.bytes_to_float3(self.tcp_recv(context, 4 * 3))

        #ack
        #tcp::recv_data((char*)&ack, sizeof(ack));
        self.tcp_send(context, bmessage_type) 

        bbox = self.get_bbox(context)
        bbox.location = ((bbbox_min[0] + bbbox_max[0]) / 2.0, (bbbox_min[1] + bbbox_max[1]) / 2.0, (bbbox_min[2] + bbbox_max[2]) / 2.0)
        #scale = max(((bbbox_max[0] - bbbox_min[0]) / (bbox.empty_display_size * 2), (bbbox_max[1] - bbbox_min[1]) / (bbox.empty_display_size * 2), (bbbox_max[2] - bbbox_min[2]) / (bbox.empty_display_size * 2)))
        scale = max(((bbbox_max[0] - bbbox_min[0]) / 2.0, (bbbox_max[1] - bbbox_min[1]) / 2.0, (bbbox_max[2] - bbbox_min[2]) / 2.0))
        bbox.scale = (scale, scale, scale)

        # def get_bbox_size(self, context):
        #     bbox = self.get_bbox(context)
        #     return bbox.empty_display_size * 2
        
        # def get_bbox_dim(self, context):
        #     bbox = self.get_bbox(context)

        #     return bbox.scale * bbox.empty_display_size * 2

        # def get_bbox_location(self, context):
        #     bbox = self.get_bbox(context)

        #     return bbox.location

        return     

    def extract_raw_part(self, context, file_data):
        import struct
        from io import BytesIO

        class ParticleData:
            def __init__(self):
                self.name = ""
                self.num_comp = 0
                self.values = []

        class RawParticles:
            def __init__(self):
                self.data = []

            def deserialize_from_vector(self, bin_data):
                stream = BytesIO(bin_data)

                # Read the size of the data vector
                data_size = struct.unpack('Q', stream.read(8))[0]  # Assuming size_t is 8 bytes (platform dependent)

                self.data = [ParticleData() for _ in range(data_size)]

                # Read each ParticleData
                for particle in self.data:
                    # Read the size of the name string
                    name_size = struct.unpack('Q', stream.read(8))[0]

                    # Read the characters of the name string
                    particle.name = stream.read(name_size).decode('utf-8')

                    # Read num_comp
                    particle.num_comp = struct.unpack('i', stream.read(4))[0]  # Assuming int is 4 bytes

                    # Read the size of the values vector
                    values_size = struct.unpack('Q', stream.read(8))[0]

                    # Read the values
                    particle.values = list(struct.unpack(f'{values_size}f', stream.read(values_size * 4)))

            def create_gnodes_instance(self, context, obj, instance_object):
                # pass
                context.view_layer.objects.active = obj
                bpy.ops.node.new_geometry_nodes_modifier()

                obj.modifiers[0].name = instance_object.name
                node_group = obj.modifiers[0].node_group
                node_group.name = instance_object.name

                # bpy.ops.node.add_node(type="GeometryNodeObjectInfo", use_transform=True)
                node_object_info = node_group.nodes.new('GeometryNodeObjectInfo')
                node_object_info.location = (-400, 600)
                node_object_info.inputs[0].default_value = instance_object
                node_object_info.inputs[1].default_value = True
                
                # bpy.ops.node.add_node(type="GeometryNodeInstanceOnPoints", use_transform=True)
                node_instance = node_group.nodes.new('GeometryNodeInstanceOnPoints')
                node_instance.location = (-50, 100)

                for p in range(len(self.data)):
                    particle = self.data[p]

                    if particle.name == "position": # skip
                        continue                

                    # bpy.ops.node.add_node(type="GeometryNodeInputNamedAttribute", use_transform=True)
                    node_in_val = node_group.nodes.new('GeometryNodeInputNamedAttribute')
                                        
                    if particle.num_comp == 1:
                        node_in_val.data_type = 'FLOAT'
                    elif particle.num_comp == 3:
                        node_in_val.data_type = 'FLOAT_VECTOR'

                    node_in_val.inputs[0].default_value = particle.name
                    node_in_val.location = (-400, 200 * (p + 1))

                node_in = node_group.nodes['Group Input']
                node_out = node_group.nodes['Group Output']

                node_group.links.clear()

                node_group.links.new(
                    node_in.outputs['Geometry'], node_instance.inputs['Points'])
                # node_group.links.new(
                #     node_in_val.outputs['Attribute'], node_instance.inputs['Rotation'])

                node_group.links.new(
                    node_object_info.outputs['Geometry'], node_instance.inputs['Instance'])

                node_group.links.new(
                    node_instance.outputs['Instances'], node_out.inputs['Geometry'])          

            def add_verts_to_mesh(self, context, empty_obj, count):
                context.view_layer.objects.active = empty_obj

                mesh = empty_obj.data
                mesh.vertices.add(count)

                for particle in self.data:
                    if particle.name == "position": # skip
                        continue

                    # mesh.attributes.new(name="loc", type='FLOAT_VECTOR', domain='POINT')
                    #mesh.attributes.new(name="rot", type='FLOAT_VECTOR', domain='POINT')
                    
                    if particle.num_comp == 1:
                        mesh.attributes.new(name=particle.name, type='FLOAT', domain='POINT')
                    elif particle.num_comp == 3:
                        mesh.attributes.new(name=particle.name, type='FLOAT_VECTOR', domain='POINT')

                mesh.update()
                mesh.validate()                      

            def create_empty_bobj(self, context, obj_collection, name):
                empty_name = '%s_empty' % name
                mesh = bpy.data.meshes.new(name=empty_name)

                mesh.update()
                mesh.validate()

                obj = bpy.data.objects.new(empty_name, mesh)
                obj_collection.objects.link(obj)

                return obj

            def create_collection(self, context):
                # Collection
                bcollection_main = bpy.data.collections.new('SPACE')
                context.scene.collection.children.link(bcollection_main)          

                return bcollection_main
            
            def extract_GN(self, context):
                self.deserialize_from_vector(file_data)    

                collection = self.create_collection(context)
                empty_obj = self.create_empty_bobj(context, collection, "SPACE")
                count = int(len(self.data[0].values) / self.data[0].num_comp)
                self.add_verts_to_mesh(context, empty_obj, count)
                
                #position
                for id in range(count):
                    for particle in self.data:
                        if particle.name == "position":
                            for j in range(3):
                                empty_obj.data.vertices[id].co[j] = particle.values[id * 3 + j]
                        else:
                            if particle.num_comp == 1:
                                empty_obj.data.attributes[particle.name].data[id].value = particle.values[id]
                            elif particle.num_comp == 3:
                                for j in range(3):
                                    empty_obj.data.attributes[particle.name].data[id].vector[j] = particle.values[id * 3 + j]

                # Add a UV sphere to the scene
                bpy.ops.mesh.primitive_uv_sphere_add(
                    radius=0.1,             # Radius of the sphere
                    location=(0, 0, 0),     # Location of the sphere
                    segments=32,            # Number of segments (horizontal slices)
                    ring_count=16           # Number of rings (vertical slices)
                )
                self.create_gnodes_instance(context, empty_obj, context.object)

        rawParticles = RawParticles()
        rawParticles.extract_GN(context)

    def extract_data(self, context):
        if context.scene.bspace_list_types_index < 0 or context.scene.bspace_list_types_index > len(context.scene.bspace_list_types):
            raise Exception("context.scene.bspace_list_types_index < 0 or context.scene.bspace_list_types_index > len(context.scene.bspace_list_types)")

        # if context.scene.view_pg_bspace.use_view_bbox == True:
        #     #    self.find_bbox(context)
        #     #else:
        #     bbox_min,bbox_max = get_view_bounds_3d(context)
        #     grid_dim = context.scene.view_pg_bspace.grid_dim
        # else:
        #     # bbox_location = context.scene.view_pg_bspace.box_object_select.location
        #     # bbox_dims = context.scene.view_pg_bspace.box_object_select.dimensions
        #     # bbox_min = (bbox_location[0] - bbox_dims[0] / 2.0, bbox_location[1] - bbox_dims[1] / 2.0, bbox_location[2] - bbox_dims[2] / 2.0)
        #     # bbox_max = (bbox_location[0] + bbox_dims[0] / 2.0, bbox_location[1] + bbox_dims[1] / 2.0, bbox_location[2] + bbox_dims[2] / 2.0)

        bbox_location = self.get_bbox_location(context)
        bbox_dims = self.get_bbox_dim(context)
        bbox_size = context.scene.view_pg_bspace.bbox_size #self.get_bbox_size(context)

        bbox_min = (bbox_location[0] - bbox_dims[0] / 2.0, bbox_location[1] - bbox_dims[1] / 2.0, bbox_location[2] - bbox_dims[2] / 2.0)
        bbox_max = (bbox_location[0] + bbox_dims[0] / 2.0, bbox_location[1] + bbox_dims[1] / 2.0, bbox_location[2] + bbox_dims[2] / 2.0)       
        
        grid_dim = context.scene.view_pg_bspace.grid_dim               

        message_type = 2
        bmessage_type = self.int_to_bytes(message_type)
        self.tcp_send(context, bmessage_type)

        context.scene.view_pg_bspace.anim_task_counter = context.scene.view_pg_bspace.anim_task_counter + 1
        #################send
        #tcp::recv_data((char*)&bbox_min[0], sizeof(float) * 3);
        bbbox_min = self.float3_to_bytes(bbox_min)
        self.tcp_send(context, bbbox_min)        

        #tcp::recv_data((char*)&bbox_max[0], sizeof(float) * 3);
        bbbox_max = self.float3_to_bytes(bbox_max)
        self.tcp_send(context, bbbox_max)           

        # tcp::recv_data((char*)&bbox_dims[0], sizeof(int) * 3);
        bvdb_size = self.int_to_bytes(grid_dim)        
        self.tcp_send(context, bvdb_size)

        # tcp::recv_data((char*)&grid_transform, sizeof(float));
        vdb_transform = 1.0 #float(context.scene.view_pg_bspace.bbox_size) / float(grid_dim)
        bvdb_transform = self.float_to_bytes(vdb_transform)        
        self.tcp_send(context, bvdb_transform)        

        # tcp::recv_data((char*)&particle_type, sizeof(int));
        particle_type = context.scene.bspace_list_types[context.scene.bspace_list_types_index].Type_id
        bparticle_type = self.int_to_bytes(int(particle_type))        
        self.tcp_send(context, bparticle_type)          

        # tcp::recv_data((char*)&block_name_id, sizeof(int));
        block_name = context.scene.bspace_list_types[context.scene.bspace_list_types_index].Block_id
        bblock_name_id = self.int_to_bytes(int(block_name))        
        self.tcp_send(context, bblock_name_id)

		# tcp::recv_data((char*)&value_convert, sizeof(int));              
        # bvalue_convert = self.int_to_bytes(int(context.scene.view_pg_bspace.value_convert))        
        # self.tcp_send(context, bvalue_convert)

        # tcp::recv_data((char*)&extracted_type, sizeof(int)); 
        bextracted_type = self.int_to_bytes(int(context.scene.view_pg_bspace.extracted_type))        
        self.tcp_send(context, bextracted_type)

		# tcp::recv_data((char*)&dense_type, sizeof(int));              
        bdense_type = self.int_to_bytes(int(context.scene.view_pg_bspace.dense_type))        
        self.tcp_send(context, bdense_type)

        bdense_norm = self.int_to_bytes(int(context.scene.view_pg_bspace.dense_norm))        
        self.tcp_send(context, bdense_norm)

        # tcp::recv_data((char*)&bbox_size, sizeof(float));
        bobject_size= self.float_to_bytes(float(bbox_size))
        self.tcp_send(context, bobject_size)

        #tcp::recv_data((char*)&particle_fix_size, sizeof(int));        
        bparticle_fix_size= self.float_to_bytes(float(context.scene.view_pg_bspace.particle_fix_size))        
        self.tcp_send(context, bparticle_fix_size)

        #tcp::recv_data((char*)&filter_min, sizeof(float));        
        bfilter_min= self.float_to_bytes(float(context.scene.view_pg_bspace.filter_min))        
        self.tcp_send(context, bfilter_min)

        #tcp::recv_data((char*)&filter_max, sizeof(float));
        bfilter_max= self.float_to_bytes(float(context.scene.view_pg_bspace.filter_max))        
        self.tcp_send(context, bfilter_max)

        bframe = self.int_to_bytes(int(context.scene.view_pg_bspace.anim_frame))        
        self.tcp_send(context, bframe)

        banim_type = self.int_to_bytes(int(context.scene.view_pg_bspace.anim_type))
        self.tcp_send(context, banim_type)

        banim_task_counter = self.int_to_bytes(int(context.scene.view_pg_bspace.anim_task_counter))
        self.tcp_send(context, banim_task_counter)
        #################recv
        file_type_id = self.bytes_to_int(self.tcp_recv(context, 4))
        #print("file_type_id:", file_type_id)
        
        file_type = file_type_items[file_type_id][1]
        #tcp::send_data((char*)&size, sizeof(size));
        file_size = self.bytes_to_int64(self.tcp_recv(context, 8))
        #tcp::send_data(file_content.data(), size);
        if file_size > 0: 
            file_data = self.tcp_recv(context, file_size)
        else:
            file_data = None

		#vdb info
        #tcp::send_data((char*)&min_value, sizeof(min_value));
        min_value = self.bytes_to_float(self.tcp_recv(context, 4))
        #tcp::send_data((char*)&max_value, sizeof(max_value));
        max_value = self.bytes_to_float(self.tcp_recv(context, 4))

        #tcp::send_data((char*)&min_value, sizeof(min_value));
        min_value_reduced = self.bytes_to_float(self.tcp_recv(context, 4))
        #tcp::send_data((char*)&max_value, sizeof(max_value));
        max_value_reduced = self.bytes_to_float(self.tcp_recv(context, 4))

        #anim
        #tcpConnection.send_data_data((char*)&frames, sizeof(frames));
        frames = self.bytes_to_int(self.tcp_recv(context, 4))

        #ack
        #tcp::recv_data((char*)&ack, sizeof(ack));
        self.tcp_send(context, bmessage_type)   

        ##########################processing
        if file_type == "RAW_PART":
            self.extract_raw_part(context, file_data)
            return

        from datetime import datetime
        dt = datetime.now().isoformat('-').replace(':', '').replace('.', '')
        fvdb = "out_%s" % dt[0:19]

        type_name = context.scene.bspace_list_types[context.scene.bspace_list_types_index].Type
        block_name = context.scene.bspace_list_types[context.scene.bspace_list_types_index].Block
        fvdb = type_name + " - " + block_name + " - " + fvdb
        fvdb = fvdb.replace(" ", "")

        pref = bspace_pref.preferences()

        if not file_data is None:
            filename = pref.local_temp_dir_path + "/" + fvdb + ".vdb"
            if file_type == "OPENVDB":
                with open(filename, 'wb') as f:
                    f.write(file_data)

            elif file_type == "NANOVDB":
                filename_nvdb = pref.local_temp_dir_path + "/" + fvdb + ".nvdb"
                with open(filename_nvdb, 'wb') as f:
                    f.write(file_data)
                
                if len(pref.nvdb_converter_path) > 0:
                    cmd = [
                        pref.nvdb_converter_path,
                        filename_nvdb,
                        filename
                    ]

                    import subprocess
                    process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                    stdout, stderr = process.communicate()

                    if process.returncode != 0:
                        if stdout:
                            print(str(stdout.decode()))
                        if stderr:
                            print(str(stderr.decode()))        

                        raise Exception("nanovdb command failed: %s" % cmd)
            elif file_type == "PATH":
                filename = file_data.decode()
                if context.scene.view_pg_bspace.replace_path_enabled:
                    filename = filename.replace(context.scene.view_pg_bspace.replace_path_orig, context.scene.view_pg_bspace.replace_path_new)

            # hide all vdbs
            try:
                for ob in bpy.data.objects:
                    if 'BSPACE' in ob:
                        ob.hide_render = True
                        ob.hide_set(True)
            except:
                # ignore
                pass

            # TODO: handle particle extraction
            if context.scene.view_pg_bspace.extracted_type == '2': # PARTICLE
                #filename = file_data.decode()
                # Handle particle extraction
                #for f in range(context.scene.view_pg_bspace.anim_start, context.scene.view_pg_bspace.anim_end + 1):
                ##self.extract_raw_part(context, file_data) # filename.replace("00000", f"{f:05d}").decode()
                raise Exception("Particle extraction not implemented yet")
            else:
                # import vdb
                bpy.ops.object.volume_import(filepath=filename, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
        else:
            if self.get_bbox_name() in bpy.data.objects:
                obj_cube = bpy.data.objects[self.get_bbox_name()]

                if obj_cube.type == 'EMPTY':
                    # If the object is linked to a collection, unlink it first
                    for collection in obj_cube.users_collection:
                        collection.objects.unlink(obj_cube)

                    # Delete the object
                    bpy.data.objects.remove(obj_cube)

                    #################### Create a new empty object
                    #bpy.ops.object.empty_add(type='CUBE')
                    # Define the vertices for the cube            
                    v = context.scene.view_pg_bspace.bbox_size
                    vertices = [
                        (-v/2.0, -v/2.0, -v/2.0),
                        (v/2.0, -v/2.0, -v/2.0),
                        (v/2.0, v/2.0, -v/2.0),
                        (-v/2.0, v/2.0, -v/2.0),
                        (-v/2.0, -v/2.0, v/2.0),
                        (v/2.0, -v/2.0, v/2.0),
                        (v/2.0, v/2.0, v/2.0),
                        (-v/2.0, v/2.0, v/2.0)
                    ]

                    # Define the edges for the cube
                    #edges = []
                    edges = [
                        (0, 1),
                        (1, 2),
                        (2, 3),
                        (3, 0),
                        (4, 5),
                        (5, 6),
                        (6, 7),
                        (7, 4),
                        (0, 4),
                        (1, 5),
                        (2, 6),
                        (3, 7)
                    ]
                    faces = []

                    # Create the new mesh and object
                    mesh = bpy.data.meshes.new(self.get_bbox_name())
                    mesh.from_pydata(vertices, edges, faces)
                    mesh.update()

                    # if hasattr(bpy.types.Scene, 'haystack_scene'):
                    #     context.scene.haystack_scene.create_bbox(context)

                    obj_cube = bpy.data.objects.new(self.get_bbox_name(), mesh)

                    bbox_size = context.scene.view_pg_bspace.bbox_size
                    bbox_size_half = bbox_size / 2.0

                    #obj_cube.empty_display_size = bbox_size_half
                    obj_cube.scale=(bbox_size_half, bbox_size_half, bbox_size_half)
                    obj_cube.location=(bbox_size_half, bbox_size_half, bbox_size_half)                    

                    # Add the object into the scene
                    #scene = bpy.context.scene
                    #scene.collection.objects.link(obj_cube)
                    # Get the active collection
                    active_collection = bpy.context.view_layer.active_layer_collection.collection

                    # Add the object into the active collection
                    active_collection.objects.link(obj_cube)

                    # # Set the origin to the center of the cube
                    # v = context.scene.view_pg_bspace.bbox_size
                    # origin_location = (v/2, v/2, v/2)

                    # # Store the current cursor location
                    # saved_cursor_location = context.scene.cursor.location.copy()

                    # # Move cursor to the desired origin position
                    # context.scene.cursor.location = origin_location

                    # # Set origin to cursor position
                    # context.view_layer.objects.active = obj_cube
                    # bpy.ops.object.origin_set(type='ORIGIN_CURSOR')

                    # # Restore cursor position
                    # context.scene.cursor.location = saved_cursor_location                    

                    obj_cube['MIN_VALUE'] = 0.0
                    obj_cube['MAX_VALUE'] = 1.0
                    self.set_vdb_shader(context, obj_cube)

                    try:
                        context.scene.cyclesphi.server_settings.mat_volume = obj_cube.data.materials[0]
                    except:
                        pass     

                # Optionally set the object as active and select it
                bpy.context.view_layer.objects.active = obj_cube
                #obj_cube.select_set(True)                

        #set view
        # if grid_dim == context.scene.view_pg_bspace.grid_dim:
        #     context.scene.bspace.set_view(context)
        #     bpy.ops.view3d.view_selected(use_all_regions=False)

        obj_vdb_new = bpy.context.view_layer.objects.active

        if not file_data is None:
            for i in range(3):
                obj_vdb_new.lock_location[i] = True
                obj_vdb_new.lock_rotation[i] = True
                obj_vdb_new.lock_scale[i] = True
        # obj_vdb_new.scale = context.scene.view_pg_bspace.box_object_select.scale

        # bbox_location = context.scene.view_pg_bspace.box_object_select.location
        # bbox_dims = context.scene.view_pg_bspace.box_object_select.dimensions
        # bbox_min = (bbox_location[0] - bbox_dims[0] / 2.0, bbox_location[1] - bbox_dims[1] / 2.0, bbox_location[2] - bbox_dims[2] / 2.0)
        # obj_vdb_new.location = bbox_min
                
        #obj_vdb_new.matrix_world = context.scene.view_pg_bspace.box_object_select.matrix_world
        #obj_vdb_new.select_set(False)

        if not file_data is None:
            obj_vdb_new.data.render.clipping = 0
            obj_vdb_new.data.render.precision = 'FULL'

            if frames > 1 and context.scene.view_pg_bspace.anim_type == '1': # ALL:
                obj_vdb_new.data.is_sequence = True
                obj_vdb_new.data.frame_duration = frames
                obj_vdb_new.data.frame_start = 1
                obj_vdb_new.data.frame_offset = -1
                obj_vdb_new.data.sequence_mode = 'REPEAT'

        # set properties
        obj_vdb_new['BSPACE'] = True

        obj_vdb_new['MIN_VALUE'] = min_value
        obj_vdb_new['MAX_VALUE'] = max_value
        obj_vdb_new['MIN_VALUE_REDUCED'] = min_value_reduced
        obj_vdb_new['MAX_VALUE_REDUCED'] = max_value_reduced

        obj_vdb_new['GRID_DIM'] = grid_dim
        #obj_vdb_new['VDB_TRANSFORM'] = context.scene.view_pg_bspace.vdb_transform
        obj_vdb_new['PARTICLE_TYPE'] = particle_type
        obj_vdb_new['BLOCK_NAME'] = block_name
        #obj_vdb_new['VALUE_CONVERT'] = context.scene.view_pg_bspace.value_convert
        obj_vdb_new['EXTRACTED_TYPE'] = context.scene.view_pg_bspace.extracted_type
        obj_vdb_new['DENSE_TYPE'] = context.scene.view_pg_bspace.dense_type
        obj_vdb_new['DENSE_NORM'] = context.scene.view_pg_bspace.dense_norm

        obj_vdb_new['BBOX_MIN'] = bbox_min
        obj_vdb_new['BBOX_MAX'] = bbox_max
        #obj_vdb_new['OBJECT_SIZE'] = context.scene.view_pg_bspace.bbox_size
        obj_vdb_new['PARTICLE_SIZE'] = context.scene.view_pg_bspace.particle_fix_size
        obj_vdb_new['FILTER_MIN'] = context.scene.view_pg_bspace.filter_min
        obj_vdb_new['FILTER_MAX'] = context.scene.view_pg_bspace.filter_max

        type_name = context.scene.bspace_list_types[context.scene.bspace_list_types_index].Type
        block_name = context.scene.bspace_list_types[context.scene.bspace_list_types_index].Block
        name = type_name + " - " + block_name
        obj_vdb_new['NAME'] = name

        #if file_data is None and hasattr(bpy.types.Scene, 'haystack_scene'):
        if not file_data is None:
            obj_vdb_new.name = name + " - " + obj_vdb_new.name

        # if context.scene.view_pg_bspace.use_view_bbox == False:
        #     context.scene.view_pg_bspace.box_object_select.select_set(False)
        #bpy.context.view_layer.objects.active = context.scene.view_pg_bspace.box_object_select

        #set new shader
        if not file_data is None or len(obj_vdb_new.data.materials) == 0:
            self.set_vdb_shader(context, obj_vdb_new)

        mat_node_tree = obj_vdb_new.data.materials[0].node_tree

        if "Map Range" in mat_node_tree.nodes:
            mat_node_tree.nodes["Map Range"].inputs['From Min'].default_value = obj_vdb_new['MIN_VALUE_REDUCED']
            mat_node_tree.nodes["Map Range"].inputs['From Max'].default_value = obj_vdb_new['MAX_VALUE_REDUCED']

        if "Math" in mat_node_tree.nodes:
            mat_node_tree.nodes["Math"].inputs[1].default_value = context.scene.view_pg_bspace.density

        # add to list
        if not file_data is None:
            item = context.scene.bspace_list_data.add()
            item.Id = len(context.scene.bspace_list_data)

            #item.Name = str(particle_type_items[int(obj_vdb_new['PARTICLE_TYPE'])][1]) + str(" - ") + str(block_name_items[int(obj_vdb_new['BLOCK_NAME'])][1])
            item.Name = obj_vdb_new.name
            item.Obj = obj_vdb_new.name

            if context.scene.bspace_list_data_index > len(context.scene.bspace_list_data) - 1:
                context.scene.bspace_list_data_index = len(context.scene.bspace_list_data) - 1

    def merge_vdb(self, context):
        # filename_nvdb = pref.local_temp_dir_path + "/" + fvdb + ".nvdb"
        # with open(filename_nvdb, 'wb') as f:
        #     f.write(file_data)
        
        pref = bspace_pref.preferences()
        if len(pref.vdb_merger_path) > 0:
            cmd = []
            cmd.append(pref.vdb_merger_path)

            # hide all vdbs
            for ob in bpy.data.objects:
                if 'BSPACE' in ob:
                    if ob.hide_get() == False:
                        filepath = bpy.path.abspath(ob.data.filepath)
                        cmd.append(filepath)
            
            out_file = cmd[-1].replace(".vdb", "_merged.vdb")
            cmd.append('-o')
            cmd.append(out_file)

            print("merging: %s" % cmd)          

            import subprocess
            process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stdout, stderr = process.communicate()

            if process.returncode != 0:
                if stdout:
                    print(str(stdout.decode()))
                if stderr:
                    print(str(stderr.decode()))        

                raise Exception("vdb command failed: %s" % cmd)
            
            # hide all vdbs
            try:
                for ob in bpy.data.objects:
                    if 'BSPACE' in ob:
                        ob.hide_render = True
                        ob.hide_set(True)
            except:
                # ignore
                pass
        
            # import vdb
            bpy.ops.object.volume_import(filepath=out_file, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))

    def vdb2histo(self, context):        
        pref = bspace_pref.preferences()
        if len(pref.vdb2histo_path) > 0:
            cmd = []
            cmd.append(pref.vdb2histo_path)

            for ob in bpy.data.objects:
                if 'BSPACE' in ob:
                    if ob.hide_get() == False:
                        filepath = bpy.path.abspath(ob.data.filepath)
                        cmd.append(filepath)
                        break
            
            out_file_png = cmd[-1].replace(".vdb", ".png")
            cmd.append('-o')
            cmd.append(out_file_png)

            cmd.append('-g')
            cmd.append('density')

            print("vdb2histo: %s" % cmd)

            import subprocess
            process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stdout, stderr = process.communicate()

            if process.returncode != 0:
                if stdout:
                    print(str(stdout.decode()))
                if stderr:
                    print(str(stderr.decode()))        

                raise Exception("vdb command failed: %s" % cmd)            
        
            # Import the PNG file
            image = bpy.data.images.load(out_file_png, check_existing=True)
            image.reload()  # Make sure it's up to date

    def colorramp_to_array(self, context, ramp, size):
        """
        Convert a color ramp to arrays of colors and alpha values.

        Parameters:
        ramp : A color ramp object that has an `evaluate` method.
        ramp_color : List to store the color values (as tuples of three floats).
        ramp_alpha : List to store the alpha values (as floats).
        size : The number of divisions in the ramp.
        """
        full_size = size + 1
        
        # Resize the output lists.
        ramp_color = []
        ramp_alpha = []

        for i in range(full_size):
            color = ramp.evaluate(float(i) / float(size)) # Evaluate the ramp.
            
            # Clamp color to 0-255 range and convert to integer
            ramp_color.append(int(max(0, min(255, color[0] * 255.0))))  # Extract and clamp R
            ramp_color.append(int(max(0, min(255, color[1] * 255.0))))  # Extract and clamp G
            ramp_color.append(int(max(0, min(255, color[2] * 255.0))))  # Extract and clamp B

            # Clamp alpha to 0-255 range and convert to integer
            ramp_alpha.append(int(max(0, min(255, color[3] * 255.0))))  # Extract and clamp Alpha

        return ramp_color, ramp_alpha            

    def vdb2png(self, context):
        pref = bspace_pref.preferences()
        if len(pref.vdb2png_path) > 0:
            cmd = []
            cmd.append(pref.vdb2png_path)

            node_tree = None

            for ob in bpy.data.objects:
                if 'BSPACE' in ob:
                    if ob.hide_get() == False:
                        filepath = bpy.path.abspath(ob.data.filepath)
                        cmd.append(filepath)
                        
                        try:
                            node_tree = ob.data.materials[0].node_tree
                        except:
                            pass

                        break

            slice_axis = context.scene.view_pg_bspace.slice_axis
            slice_nr = context.scene.view_pg_bspace.slice_nr

            out_file_png = cmd[-1].replace(".vdb", "_" + str(slice_axis) + "_" + f"{slice_nr:05d}" + ".png")
            cmd.append('-o')
            cmd.append(out_file_png)

            cmd.append('-g')
            cmd.append('density')

            cmd.append('--axis')
            cmd.append(str(slice_axis))

            cmd.append('--slice')
            cmd.append(str(slice_nr))

            node = None
            if node_tree is not None and "Color Ramp" in node_tree.nodes:
                node = node_tree.nodes["Color Ramp"]
            elif node_tree is not None and "ColorRamp" in node_tree.nodes:
                node = node_tree.nodes["ColorRamp"]

            if node is not None and node.color_ramp is not None:
                # Convert the color ramp to arrays
                ramp_values, ramp_alpha = self.colorramp_to_array(context, node.color_ramp, 256)
                ramp_values_str = ','.join(map(str, ramp_values))
                cmd.append('--ramp')
                cmd.append(ramp_values_str)

            print("vdb2png: %s" % cmd)

            import subprocess
            process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stdout, stderr = process.communicate()

            if process.returncode != 0:
                if stdout:
                    print(str(stdout.decode()))
                if stderr:
                    print(str(stderr.decode()))        

                raise Exception("vdb command failed: %s" % cmd)            
        
            # Import the PNG file
            image = bpy.data.images.load(out_file_png, check_existing=True)
            #image.source = 'SEQUENCE'  # Enable sequence mode
            image.reload()  # Make sure it's up to date  

        pass

    def recvall(self, channel, expected_length):
        data = bytearray()
        while len(data) < expected_length:
            chunk = channel.recv(expected_length - len(data))
            if not chunk:
                break  # Connection closed
            data.extend(chunk)
        return data               

    def tcp_send(self, context, data):
        if self.connection_channel is None:
            return
                
        self.connection_channel.sendall(data)

    def tcp_recv(self, context, size):
        if self.connection_channel is None:
            return None
                
        return self.recvall(self.connection_channel, size)

    def connect(self, context):
        bspace_remote.connect_to_server(context)

    def disconnect(self, context):
        if self.connection_channel is None:
            return
        
        #pref = bspace_pref.preferences()
        #if pref.bspace_enable_ssh == True:
        self.connection_channel.close()
            # self.connection_ssh2.close()
            # self.connection_ssh1.close()

        #self.connection_forward_server.shutdown()
        #self.connection_thread.join()

        # self.connection_ssh1 = None
        # self.connection_ssh2 = None
        #self.connection_thread = None
        #self.connection_forward_server = None    
        self.connection_channel = None

class BSpacePanel:
    """Creates a BSPACE Panel"""
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "BSPACE"

class BSPACE_PT_PANEL(BSpacePanel, bpy.types.Panel):
    """Creates a BSPACE Panel"""
    bl_label = "BSPACE"
    bl_idname = "BSPACE_PT_PANEL"

    def draw(self, context):
        # layout = self.layout
        # scene = context.scene
        # pref = bspace_pref.preferences()
        pass

class BSPACE_PT_local_settings(BSpacePanel, bpy.types.Panel):
    bl_label = "Local Settings"
    bl_parent_id = "BSPACE_PT_PANEL"
    #bl_idname = "BSPACE_PT_local_settings"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        pref = bspace_pref.preferences()

        col = layout.column()
        # col.prop(scene.view_pg_bspace, "local_data_dir_path")
        col.prop(pref, "local_temp_dir_path")
        #col.operator("bspace.load_data")  

class BSPACE_PT_remote_settings(BSpacePanel, bpy.types.Panel):
    bl_label = "Remote Settings"
    bl_parent_id = "BSPACE_PT_PANEL"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout
        pref = bspace_pref.preferences()

        col = layout.column()
        col.prop(pref, "bspace_server_name", text="Server")
        col.prop(pref, "bspace_port", text="Port")   

class BSPACE_PT_environment(BSpacePanel, bpy.types.Panel):
    bl_label = "Environment"
    bl_parent_id = "BSPACE_PT_PANEL"
    #bl_idname = "BSPACE_PT_environment"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout

        col = layout.column()
        col.operator("bspace.set_world_environment")

class BSPACE_PT_geometry(BSpacePanel, bpy.types.Panel):
    bl_label = "Geometry"
    bl_parent_id = "BSPACE_PT_PANEL"
    #bl_idname = "BSPACE_PT_extract"
    bl_options = {'DEFAULT_CLOSED'}

    # @classmethod
    # def poll(cls, context):
    #     return context.scene.bspace.connection_channel

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        pref = bspace_pref.preferences()   

        box = layout.box()
        box.label(text="Geometry")
        col = box.column()            
        col.prop(scene.view_pg_bspace, "bbox_size")
        col.prop(scene.view_pg_bspace, "grid_dim")            
        #col.prop(scene.view_pg_bspace, "vdb_transform")
        #col = box.column()            
        #col.prop(scene.view_pg_bspace, "use_view_bbox") 
        #col.operator("bspace.find_types")

        # if scene.view_pg_bspace.use_view_bbox == False:
        #     col = box.column()
        #     # col.enabled = False
        #     # col.prop(scene.view_pg_bspace, "box_object_main") 
        #     # col = box.column()
        #     # col.enabled = False            
        #     # col.prop(scene.view_pg_bspace, "box_object_select")
        #     # col = box.column()
        col.operator("bspace.create_bbox")

class BSPACE_PT_connection(BSpacePanel, bpy.types.Panel):
    bl_label = "Connection"
    bl_parent_id = "BSPACE_PT_PANEL"
    #bl_idname = "BSPACE_PT_connection"
    #bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        pref = bspace_pref.preferences()

        col = layout.column()                       

        if scene.bspace.connection_channel is None:
            col.operator("bspace.connect")
        else:
            col.operator("bspace.disconnect")     

class BSPACE_PT_extract(BSpacePanel, bpy.types.Panel):
    bl_label = "Extract"
    bl_parent_id = "BSPACE_PT_PANEL"
    #bl_idname = "BSPACE_PT_extract"
    #bl_options = {'DEFAULT_CLOSED'}

    # @classmethod
    # def poll(cls, context):
    #     return context.scene.bspace.connection_channel

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        #pref = bspace_pref.preferences()   

        box = layout.box()
        box.label(text="Extract")
        # col = box.column()
        # col.enabled = False
        # col.prop(scene.view_pg_bspace, "snapshot_size")
        col = box.column()
        #col.prop(scene.view_pg_bspace, "particle_type")
        #col.prop(scene.view_pg_bspace, "block_name")

        #table
        col.template_list("BSPACE_UL_ExtractedTypes", "", scene, "bspace_list_types", scene, "bspace_list_types_index")
        #col.operator('bspace.select_from_list').index = context.scene.bspace_list_data_index          

        #col.prop(scene.view_pg_bspace, "value_convert")
        col.prop(scene.view_pg_bspace, "extracted_type")
        col.prop(scene.view_pg_bspace, "dense_type")
        col.prop(scene.view_pg_bspace, "dense_norm")
        col.prop(scene.view_pg_bspace, "particle_fix_size", text="Radius factor")
        col.prop(scene.view_pg_bspace, "filter_min", text="Min. Value")
        col.prop(scene.view_pg_bspace, "filter_max", text="Max. Value")
        col.prop(scene.view_pg_bspace, "density", text="Mat. Density")
        col.prop(scene.view_pg_bspace, "anim_type", text="Anim Type")
        
        if scene.view_pg_bspace.anim_type != '0': # Anim
            col.prop(scene.view_pg_bspace, "anim_frame", text="Frame")
            col.prop(scene.view_pg_bspace, "anim_task_counter", text="Task Counter")
            col = box.column()
            col.enabled = False
            col.prop(scene.view_pg_bspace, "anim_start", text="Start")
            col.prop(scene.view_pg_bspace, "anim_end", text="End")

        if scene.view_pg_bspace.anim_type == '1': # ALL Path
            col = box.column()
            col.prop(scene.view_pg_bspace, "replace_path_enabled", text="Replace Path")
            if scene.view_pg_bspace.replace_path_enabled:
                col.prop(scene.view_pg_bspace, "replace_path_orig", text="Orig Path")
                col.prop(scene.view_pg_bspace, "replace_path_new", text="New Path")            

        box = layout.box()
        col = box.column()        
        col.operator("bspace.extract_data")
        col.prop(scene.view_pg_bspace, "register_export", text="Register")

        # box = layout.box()
        # box.label(text="Utils")
        # col = box.column()
        # col.operator("bspace.calc_bbox_viewport")

class BSPACE_PT_data(BSpacePanel, bpy.types.Panel):
    bl_label = "Data"
    bl_parent_id = "BSPACE_PT_PANEL"
    #bl_idname = "BSPACE_PT_data"
    #bl_options = {'DEFAULT_CLOSED'}

    # @classmethod
    # def poll(cls, context):
    #     return context.scene.bspace.connection_channel    

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        pref = bspace_pref.preferences()

        try:
            box = layout.box()
            box.label(text="Material") 
            col = box.column()
            obj_vdb = bpy.context.view_layer.objects.active
            col.prop(obj_vdb.data.materials[0].node_tree.nodes["Math"].inputs[1], "default_value", text="Density")  
        except:
            pass 

        #table
        box = layout.box()
        col = box.column()
        col.operator("bspace.find_bbox")
        col.operator("bspace.select_bbox")
        #col.operator("bspace.move_to_click")
        col.operator("bspace.move_to_cursor")

        box = layout.box()
        box.label(text="Select VDB")
        col = box.column()
        col.template_list("BSPACE_UL_ExtractedData", "", scene, "bspace_list_data", scene, "bspace_list_data_index")
        #col.operator('bspace.select_from_list').index = context.scene.bspace_list_data_index

        box = layout.box()
        box.label(text="Merge VDB")
        col = box.column()
        col.operator("bspace.merge_deselect_all")
        col.template_list("BSPACE_UL_MergeData", "", scene, "bspace_list_data", scene, "bspace_list_data_index")
        col.operator("bspace.merge_vdb")

        box = layout.box()
        box.label(text="Histogram")
        col = box.column()
        col.operator("bspace.vdb2histo")        

        box = layout.box()
        box.label(text="2D Image")
        col = box.column()
        col.prop(scene.view_pg_bspace, "slice_axis", text="Slice Axis")
        col.prop(scene.view_pg_bspace, "slice_nr", text="Slice Number")
        col.operator("bspace.vdb2png")

#####################################################################################################################    
addon_keymaps = []

# Registering
classes = [
    BSPACE_PG_SETTINGS,

    BSPACE_PG_ExtractedData,
    BSPACE_UL_ExtractedData,

    BSPACE_UL_MergeData,

    BSPACE_PG_ExtractedTypes,
    BSPACE_UL_ExtractedTypes,      

    BSPACE_OT_LoadData,
    BSPACE_OT_Connect,
    BSPACE_OT_Disconnect,
    BSPACE_OT_ExtractData,
    BSPACE_OT_SelectBBOX,
    BSPACE_OT_FindBBOX,
    BSPACE_OT_CreateBBox,
    BSPACE_OT_FindTypes,    
    BSPACE_OT_SetWorldEnvironment,   
    BSPACE_OT_select,
    BSPACE_OT_CalcBBOXViewport,
    BSPACE_OT_Zoom,
    BSPACE_OT_update_extracted_data,

    BSPACE_OT_merge_select_data,

    BSPACE_OT_move_to_click,
    BSPACE_OT_move_to_cursor,

    BSPACE_OT_MergeVDB,
    BSPACE_OT_MergeDeselectAll,

    BSPACE_OT_HistoVDB,
    BSPACE_OT_VDB2PNG,

    BSPACE_PT_PANEL,
    BSPACE_PT_local_settings,
    BSPACE_PT_remote_settings,    
    BSPACE_PT_environment,
    BSPACE_PT_geometry,
    BSPACE_PT_connection,
    BSPACE_PT_extract,
    BSPACE_PT_data,
]

def register():
    for cls in classes:
        bpy.utils.register_class(cls)

    bpy.types.Scene.view_pg_bspace = bpy.props.PointerProperty(type=BSPACE_PG_SETTINGS)

    bpy.types.Scene.bspace = BSPACE()

    bpy.types.Scene.bspace_list_data = bpy.props.CollectionProperty(type=BSPACE_PG_ExtractedData)
    bpy.types.Scene.bspace_list_data_index = bpy.props.IntProperty(default=-1)

    bpy.types.Scene.bspace_list_types = bpy.props.CollectionProperty(type=BSPACE_PG_ExtractedTypes)
    bpy.types.Scene.bspace_list_types_index = bpy.props.IntProperty(default=-1)

    # Add the hotkey
    wm = bpy.context.window_manager
    kc = wm.keyconfigs.addon
    if kc:
        km = wm.keyconfigs.addon.keymaps.new(name='3D View', space_type='VIEW_3D')
        kmi = km.keymap_items.new("bspace.zoom", type='W', value='PRESS', ctrl=True)
        addon_keymaps.append((km, kmi))    

def unregister():
    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)

    delattr(bpy.types.Scene, "bspace")
    delattr(bpy.types.Scene, "view_pg_bspace")

    delattr(bpy.types.Scene, "bspace_list_data")
    delattr(bpy.types.Scene, "bspace_list_data_index")

    delattr(bpy.types.Scene, "bspace_list_types")
    delattr(bpy.types.Scene, "bspace_list_types_index")

    # Remove the hotkey
    for km, kmi in addon_keymaps:
        km.keymap_items.remove(kmi)
    addon_keymaps.clear()              
