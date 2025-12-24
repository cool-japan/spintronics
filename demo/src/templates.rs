//! Askama templates for web pages

use askama::Template;

#[derive(Template)]
#[template(path = "index.html")]
pub struct IndexTemplate {}

#[derive(Template)]
#[template(path = "llg.html")]
pub struct LlgTemplate {}

#[derive(Template)]
#[template(path = "spin_pumping.html")]
pub struct SpinPumpingTemplate {}

#[derive(Template)]
#[template(path = "materials.html")]
pub struct MaterialsTemplate {}

#[derive(Template)]
#[template(path = "skyrmion.html")]
pub struct SkyrmionTemplate {}
