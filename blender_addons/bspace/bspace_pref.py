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
import bpy

ADDON_NAME = 'bspace'
#######################BSpacePreferences#########################################

class BSpacePreferences(bpy.types.AddonPreferences):
    bl_idname = ADDON_NAME

    bspace_port: bpy.props.IntProperty(
        name="Port",
        min=0,
        max=65565,
        default=5000
    ) # type: ignore

    
    bspace_server_name: bpy.props.StringProperty(
        name="Server",
        default="localhost"
    ) # type: ignore

    local_temp_dir_path: bpy.props.StringProperty(
        name="Temp", 
        subtype='DIR_PATH', 
        default=''
        )  # type: ignore

    nvdb_converter_path: bpy.props.StringProperty(
        name="NVDB Converter",
        subtype='FILE_PATH',
        default=""
    ) # type: ignore

    vdb_merger_path: bpy.props.StringProperty(
        name="VDB Merger",
        subtype='FILE_PATH',
        default=""
    ) # type: ignore

    vdb2histo_path: bpy.props.StringProperty(
        name="VDB to Histogram",
        subtype='FILE_PATH',
        default=""
    ) # type: ignore

    vdb2png_path: bpy.props.StringProperty(
        name="VDB to PNG",
        subtype='FILE_PATH',
        default=""
    ) # type: ignore     

    def draw(self, context):
        layout = self.layout

        box = layout.box()
        box.label(text='Local Settings:')
        col = box.column()
        col.prop(self, "local_temp_dir_path", text="Temp")       

        box = layout.box()
        box.label(text='BSpace TCP Server:')
        col = box.column()
        col.prop(self, "bspace_server_name", text="Server")
        col.prop(self, "bspace_port", text="Port")

        box = layout.box()
        box.label(text='Other:')
        col = box.column()
        col.prop(self, "nvdb_converter_path")
        col.prop(self, "vdb_merger_path")
        col.prop(self, "vdb2histo_path")
        col.prop(self, "vdb2png_path")

def ctx_preferences():
    try:
        return bpy.context.preferences
    except AttributeError:
        return bpy.context.user_preferences

def preferences() -> BSpacePreferences:
    return ctx_preferences().addons[ADDON_NAME].preferences

def register():
    bpy.utils.register_class(BSpacePreferences)

def unregister():
    bpy.utils.unregister_class(BSpacePreferences)