#include "ui.h"

status_t ui_init(ui_t* ui, canvas_t* canvas)
{
	status_t status = font_init(&ui->font);
	if (STATUS_OK != status)
	{
		return status;
	}

	ui->canvas = canvas;

	return STATUS_OK;
}

void ui_draw(ui_t* ui, float upscaling_factor)
{
	for (int i = 0; i < ui->text_count; ++i)
	{
		// TODO: text_t needs to be scaled. Or how do I do it not scaled.
		text_draw(ui->canvas, &ui->text[i], &ui->font, upscaling_factor);
	}
}

void ui_destroy(ui_t* ui)
{

}
